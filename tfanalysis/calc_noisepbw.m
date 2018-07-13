function NPBW = calc_noisepbw(win)
% CALC_NOISEPBW returns the equivalent noise power bandwidth of the specified
% window
%
% Ref: Haykin (2001) Communication Systems, 4th Ed., John Wiley & Sons, p. 722

% todo:
%    Add ability to numerically evaluate arbitrary window function

% force input into window structure
if ~isstruct(win)
    param = win;
    clear win;      % move input parameter so we can redefine as a struct
    if ischar(param)
        win.name = param;
        if exist(win.name,'file')
            win.fn = str2func(win.name);
        else
            warning('CALC_NOISEPBW:main','Window name not found on MATLAB path')
        end
    elseif isa(param,'function_handle')
        win.name = func2str(param);
        win.fn = param;
        if ~exist(win.name,'file')
            warning('CALC_NOISEPBW:main','Window function must be a valid function name or handle')
        end
    elseif isa(param,'double')
        win.name = 'unknown';
        win.fn = @(x) param;
    else
        error('CALC_NOISEPBW:main','Could not determine window NPBW.  Input parameter not a valid type.')
    end
    
end

% % attempt to use window function, otherwise use lookup table
% if isfield(win,'fn')
%     NPBW = calc_npbw(win.fn);            % window noise power bandwidth (calculation)
% elseif isfield(win,'name')
if exist(win.name,'file')        % temporarily use only a lookup table until calculation implemented
    NPBW = calc_npbw(win.fn);
else
    NPBW = lookup_npbw(win.name);        % window noise power bandwidth (lookup table)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NPBW = lookup_npbw(winname)
% LOOKUP_NPBW returns the equivalent noise power bandwidth based upon a
% lookup table of common windows

switch lower(winname)
    % windows available in MATLAB signal processing toolbox
    case 'rectwin'
        NPBW = 1;
    case 'triang'
        NPBW = 1.33;
    case 'hanning'
        NPBW = 1.5;
    case 'hann'
        NPBW = 1.5;
    case 'hamming'
        NPBW = 1.362826;
    case 'blackmanharris'
        NPBW = 2.004353;
    case 'flattopwin'
        NPBW = 3.822108760;
 
    % unavailable in MATLAB signal processing toolbox (Audio Precision builtin windows)
    case 'equiripple'
        NPBW = 2.631905018;
    case 'rifevincent4'
        NPBW = 2.310000000;
    case 'rifevincent5'
        NPBW = 2.626530612;
    case 'gaussian'     % alpha = ?
        NPBW = 2.215349682;
        
    otherwise
        warning('CALC_NOISEPBW:lookup','Unknown window name - Default to 1')
        NPBW = 1;
end


function B = calc_npbw(win)
% CALC_NPBW returns the equivalent noise power bandwidth for any window

N = 101;
fs = N;
nfft = 2^14;
x = win(N);
df = fs/nfft;       % spectral resolution

% compute FFT
f = (0:nfft-1)*df;
y = fft(x,nfft)./N;
yy = abs(y).^2;

% compute equivalent noise power bandwidth
n = trapz(yy)*df;
d = max(yy);

B = n/d;
