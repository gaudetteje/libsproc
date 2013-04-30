function fd = calc_spectrum(x,varargin)
%CALC_SPECTRUM   Calculates the amplitude or power spectrum of a time series signal
%
% FD = CALC_SPECTRUM(X) calculates the amplitude spectrum for a time-series
%     signal, X, with the sampling rate normalized to 1 Hz
% FD = CALC_SPECTRUM(X,FS)  specifies the sampling rate, FS
% FD = CALC_SPECTRUM(TS)  alternatively accepts the time series struct, TS
% FD = CALC_SPECTRUM(...,NFFT,WIN,NAVG,OVERLAP,BAND,ACDC,UNITS) applies the
%     optional parameters to spectral analysis
% FD = CALC_SPECTRUM(1)  returns an empty struct - useful for initialization
%     of large arrays
%
% TS is a struct with the following fields:
%   .data  LxN matrix for L samples and N channels [Volts or counts]
%   .fs    scalar sampling rate [Hz]
%
% Note:  multiple data channels, N, may be entered as an LxN matrix
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
% UNITS     'Vpk'         Specify the units of the magnitude (see below)
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
% Email:    jason.e.gaudette@navy.mil
% Date:     20130430
%

%% Handle data struct

% force data into time series structure
if isstruct(x)
    ts = x;
else
    ts.data = x;
    if nargin > 1
        ts.fs = varargin{1};
        varargin(1) = [];          % remove entry from optional input cell array
    end
end

% force data to be a column vector if only a single channel
if size(ts.data,1) == 1
    ts.data = ts.data(:);
end
if size(ts.data,2) > size(ts.data,1)
    warning('CALC_SPECTRUM:size','Number of channels exceeds number of data samples.  Verify data matrix dimensions.')
end

% verify structure contents
if ~isfield(ts,'data')
    error('Time series struct must contain "data" field')
end
if ~isfield(ts,'fs')
    ts.fs = 1;
end

%% Handle optional input parameters
%[nfft, win, navg, overlap, band, acdc, units]
optargs =  {4096, @hamming, 1, 0.5, 'half', 'dc', 'Vpk'};
nArgs = length(varargin);
optargs(1:nArgs) = varargin;                                % overwrite defaults
[nfft, win, navg, overlap, band, acdc, units] = optargs{:}; % assign parameters


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


