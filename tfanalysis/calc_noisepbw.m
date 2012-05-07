function NPBW = calc_noisepbw(win)
% CALC_NOISEPBW returns the equivalent noise power bandwidth of the specified
% window


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
        win.fn = param;
        win.name = func2str(win.fn);
        if ~exist(win.name,'file')
            warning('CALC_NOISEPBW:main','Window function must be a valid function name or handle')
        end
    else
        error('CALC_NOISEPBW:main','Could not determine window NPBW.  Input parameter not a valid type.')
    end
    
end

% % attempt to use window function, otherwise use lookup table
% if isfield(win,'fn')
%     NPBW = calc_npbw(win.fn);            % window noise power bandwidth (calculation)
% elseif isfield(win,'name')
if 1        % temporarily use only a lookup table until calculation implemented
    NPBW = lookup_npbw(win.name);        % window noise power bandwidth (lookup table)
else
    error('CALC_NOISEPBW:main','Window structure must contain fields with a function handle "fn" and/or function name "name".')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NPBW = lookup_npbw(winname)
% LOOKUP_NPBW returns the equivalent noise power bandwidth based upon a
% lookup table of common windows

switch lower(winname)
    % available in MATLAB signal processing toolbox
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


function NPBW = calc_npbw(winfn)
% CALC_NPBW returns the equivalent noise power bandwidth for any window
% callable by MATLAB

% N = 20;
% nfft = 4096;
% x = ones(N,1);
% y0 = abs(fft(x,nfft));
% y1 = abs(fft(x.*winfn(N),1024));

%TBD
NPBW = 1;