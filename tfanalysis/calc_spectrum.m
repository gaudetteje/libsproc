function fd = calc_spectrum(ts,varargin)
%CALC_SPECTRUM   Calculates the amplitude or power spectrum of a time series signal
%
% FD = CALC_SPECTRUM(TS) calculates the amplitude spectrum for a time-series
%     signal
% FD = CALC_SPECTRUM(TS,NFFT,WIN,NAVG,OVERLAP,BAND,ACDC,UNITS) applies the
%     optional parameters to spectral analysis
% FD = CALC_SPECTRUM(1)  returns an empty struct - useful for initialization
%     of large arrays
%
% TS is a struct with the following fields:
%   .data  LxN matrix for L samples and N channels [Volts or counts]
%   .fs    scalar sampling rate [Hz]
%
% TS can also be an LxN matrix, but data is normalized to a 1Hz sampling rate
%
% Optional Input Parameters:
% <name>    <default>     <description>
% NFFT      4096          FFT points
% WIN       @hamming      Any valid window function name or handle
%                           found in the MATLAB path
% NAVG      4             Number of averages to use
% OVERLAP   0.5           Percent overlap for use when averaging
% BAND      'half'        Specify a 'full' or 'half' spectrum
% ACDC      'dc'          Specify 'ac' or 'dc' coupling to remove DC offset
% UNITS     'Vrms'        Specify the units of the magnitude (see below)
%
% UNITS are initially in 'Vpk' and are then converted using CONVERT_SPECTRUM.  
%
% Output structure:
%
% FD.
%    mag        - Nx1 array of magnitude [arbitrary units]
%    magdb      - Nx1 array of magnitude [dB]
%    phase      - Nx1 array of phase information [deg]
%    freq       - Nx1 array of frequency axis [Hz]
%
%    fs         - sampling rate [Hz]
%    nfft       - number of points in FFT
%    winname    - name of window function applied
%    resolution  - frequency bin resolution
%    navg       - Number of averages
%    units      - Y-Axis Units
%    stamp      - Timestamp of creation
% 
% See also CONVERT_SPECTRUM, FFT, PWELCH

% Author:   Jason Gaudette
% Company:  Naval Undersea Warfare Center (Newport, RI)
% Phone:    401.832.6601
% Email:    gaudetteje@npt.nuwc.navy.mil
% Date:     20101230
%


% default parameters
nfft = 4096;
win = @hamming;
navg = 1;
overlap = 0.5;
band = 'half';
acdc = 'dc';
units = 'Vpk';


% parse input parameters
switch nargin
    case 8
        units = varargin{7};
        acdc = lower(varargin{6});
        band = lower(varargin{5});
        overlap = varargin{4};
        navg = varargin{3};
        win = lower(varargin{2});
        nfft = varargin{1};
    case 7
        acdc = lower(varargin{6});
        band = lower(varargin{5});
        overlap = varargin{4};
        navg = varargin{3};
        win = lower(varargin{2});
        nfft = varargin{1};
    case 6
        band = lower(varargin{5});
        overlap = varargin{4};
        navg = varargin{3};
        win = lower(varargin{2});
        nfft = varargin{1};
    case 5
        overlap = varargin{4};
        navg = varargin{3};
        win = lower(varargin{2});
        nfft = varargin{1};
    case 4
        navg = varargin{3};
        win = lower(varargin{2});
        nfft = varargin{1};
    case 3
        win = lower(varargin{2});
        nfft = varargin{1};
    case 2
        nfft = varargin{1};
end



%% Handle input parameters and data

% force into time series structure
if ~isfield(ts,'data')
    data = ts;
    clear ts;
    ts.data = data;
    clear data;
end

% force ts to be a column vector
if size(ts.data,2) > size(ts.data,1)
    ts.data = ts.data.';
end

% get both window name and function handle, when available
if ischar(win)
    winname = win;
    if exist(winname,'file')
        winfn = str2func(win);
    else
        warning('CALC_SPECTRUM:window','Window name not found on MATLAB path')
    end
elseif isa(win,'function_handle')
    winfn = win;
    winname = func2str(win);
    if ~exist(winname,'file')
        error('CALC_SPECTRUM:window','Window parameter must be a valid function name or handle')
    end
else
    error('CALC_SPECTRUM:window','Window parameter must be a valid function name or handle')
end

% if needed, determine sampling frequency from time vector
if ~isfield(ts,'fs')
    if isfield(ts,'time')
        ts.fs = 1/(ts.time(2)-ts.time(1));
    else
        warning('CALC_SPECTRUM:fs','Could not determine sampling rate; using fs=1Hz')
        ts.fs = 1;
    end
end

nChan = size(ts.data,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Process input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% condition time series data by removing DC component
if strcmp(acdc,'ac')
    ts.data = ts.data-mean(mean(ts.data));
end

% determine signal length:  case 1 (L < NFFT) or case 2 (NFFT >= L)
N = min(length(ts.data),nfft);
if navg>1 && N<nfft
    error('Not enough data samples to average')
end

% iterate FFT over blocks of data, averaging resultant PSD
Y = zeros(nfft,nChan,navg);   % preallocate 3D matrix
for n = 1:navg
    % update data index
    idx = (1:N) + floor((n-1)*N*overlap);
    
    % construct and apply window
    win = window(winfn,N);
    x = ts.data(idx,:) .* repmat(win,1,nChan);
    
    % perform FFT and normalize to number of samples (non-zero padded)
    Y(:,:,n) = fft(x, nfft)./N;
end

% calculate magnitude response and average
fd.mag = abs(Y);
fd.mag = mean(fd.mag,3);
fd.phase = unwrap(angle(Y));        % calculate and unwrap phase response

% generate reference frequency vector and truncate/shift band
switch band
    case 'full'
        fd.freq = ts.fs.*linspace(-0.5, 0.5, nfft+1).';
        fd.freq(end) = [];
        fd.mag = fftshift(fd.mag);
    case 'half'
        fd.freq = ts.fs.*linspace(0, 0.5, nfft/2+1).';
        fd.freq(end) = [];
        fd.mag = 2.*fd.mag(1:floor(nfft/2),:);   % scale spectrum for real signals
        if N==nfft
            % correct the DC level (also need to figure in window effect for N<NFFT)
            fd.mag(1,:) = fd.mag(1,:)./2;
        else
            % Need to correct DC level for Sinc width & window effect
            warning('CALC_SPECTRUM:dc','DC component not corrected (2x true DC level)')
        end
    otherwise
        error('Invalid entry for parameter "band"')
end


% calculate bin resolution, window scaling factor, and noise power bandwidth
binres = ts.fs/nfft;                    % bin resolution
SF = mean(win);                         % window scaling factor
NPBW = calc_noisepbw(winname);          % window equivalent noise power bandwidth

% amplitude spectral density with window scaling factor correction
fd.mag = fd.mag ./ SF;                  % sqrt(PSD/SF^2)

% convert to dB [dBVpk/rtbin]
fd.magdb = zeros(size(fd.mag));      % placeholder -> dB conversion performed in convert_spectrum()



%% Now add supplementary information

% static parameters
fd.fs = ts.fs;
fd.nfft = nfft;
fd.binres = binres;
fd.win.name = winname;
fd.win.fh = winfn;
fd.win.SF = SF;
fd.win.NPBW = NPBW;
fd.navg = navg;
fd.units = 'Vpk';               % default FFT units
fd.stamp = datestr(now());



%% Finally, convert spectrum to desired units
fd = convert_spectrum(fd,units);


