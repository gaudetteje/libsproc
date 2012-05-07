function x = gen_ifpulse_multi(fs,varargin)
% GEN_IFPULSE_MULTI  constructs a multi component synthetic waveform based on
% separate mono-component pulses defined by desired IF and IA
%
% Inputs:
%   fs
%   PT           - array of structs containing individual pulse definitions
%     .time      - start and stop times for each component [sec]
%     .freq      - start and stop frequencies for each component [Hz]
%     .modFn     - pulse modulation type ['cw','lfm','sfm','hfm']
%     .winFn     - amplitude shading function [window function string or handle]
%     .gain      - amplitude scaling factor for each component
%     .phase     - initial phase for each component [radians]
%
% Outputs:
%   x            - column vector of time series amplitude
%
% See also GEN_IFPULSE

% To Do:
% - add band-limited noise generation process
% - include spline functions in IF definition
% - include spline functions in windowing definition
% - allow optional parameters for windowing functions
% - create 'default' object using OOP Class - returns with basic parameters



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input parameter defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P.modFn = {'lfm'};
P.winFn = {@rectwin};       %% need to be able to pass optional window params.  Also need to accept spline 'functions'
P.gain = 1;
P.phase = 0;
fMod = 2500;     % add this value as optional 'modOpt' paramater in P struct?


% assign input parameters to working variables
assert((nargin >= 2) && (nargin <= 7), 'Incorrect number of parameters entered')
switch nargin
    case 2
        fields = fieldnames(varargin{1});
        for i = 1:length(fields)
            P.(fields{i}) = varargin{1}.(fields{i});
        end
        
    otherwise
        fields = {'time','freq','modFn','winFn','gain','phase'};
        for i = 1:nargin-1
            P.(fields{i}) = varargin{i-1};
        end

end
clear fields

assert(isfield(P,'time') && isfield(P,'freq'), ...
    'Missing fields:  P.time and P.freq must be entered in the struct')

% initialize time series signal
x = zeros(round(max(max(P.time))*fs)+1,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterate over each pulse component
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = size(P.freq,2);
for h = 1:M

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % error check input paramaters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % check type constraints
    assert(isa(P.time, 'numeric'), ...
        'P.time must be a numeric array or matrix: found "%s"\n', class(P.time))
    assert(isa(P.freq, 'numeric'), ...
        'P.freq must be a numeric array or matrix: found "%s"\n', class(P.freq))
    assert(isa(P.modFn, 'cell') || isa(P.modFn, 'char'), ...
        'P.modFn must be a string or cell array: found "%s"\n', class(P.modFn))
    assert(isa(P.winFn, 'cell') || isa(P.winFn, 'char'), ...
        'P.winFn must be a string or cell array: found "%s"\n', class(P.winFn))
    assert(isa(P.gain, 'numeric'), ...
        'P.gain must be a numeric scalar or array: found "%s"\n', class(P.gain))
    assert(isa(P.phase, 'numeric'), ...
        'P.phase must be a numeric scalar or array: found "%s"\n', class(P.phase))
    
    % check size constraints
    assert(all(size(P.time) == [2 1]) || all(size(P.time) == [2 M]), ...
        'P.time must contain a 2x1 or 2xM matrix of start and stop times: found %dx%d; M=%d\n', ...
        size(P.time,1), size(P.time,2), M)
    assert(all(size(P.freq) == [2 M]), ...
        'P.freq must contain a 2xM matrix of start and stop frequencies: found %dx%d; M=%d\n', ...
        size(P.freq,1), size(P.freq,2), M)
    assert(all(size(P.modFn) == [1 1]) || all(size(P.modFn) == [1 M]) || isa(P.modFn,'char'), ...
        'P.modFn must contain a 1x1 or 1xM cell of char arrays or a string: found %dx%d; M=%d\n', ...
        size(P.modFn,1), size(P.modFn,2), M)
    assert(all(size(P.winFn) == [1 1]) || all(size(P.winFn) == [1 M]) || isa(P.winFn,'char'), ...
        'P.winFn must contain a 1x1 or 1xM cell of char arrays, a string, or a function handle: found %dx%d; M=%d\n', ...
        size(P.winFn,1), size(P.winFn,2), M)
    assert(all(size(P.gain) == [1 1]) || all(size(P.gain) == [1 M]), ...
        'P.gain must contain a 1x1 or 1xM array of scalar gain values: found %dx%d; M=%d\n', ...
        size(P.gain,1), size(P.gain,2), M)
    assert(all(size(P.phase) == [1 1]) || all(size(P.phase) == [1 M]), ...
        'P.phase must contain a 1x1 or 1xM array of phase values: found %dx%d; M=%d\n', ...
        size(P.phase,1), size(P.phase,2), M)

    %%% add additional sanity checks here


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parse through input parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % extract specified time information
    if (size(P.time,2) > 1)
        T = P.time(:,h);
    else
        T = P.time;
    end
    t0 = T(1);
    t1 = T(2);
    T = t1-t0;
    assert(T > 0, 'P.time must increase in value')


    % extract specified frequency information
    F = P.freq(:,h);
    f0 = F(1);
    f1 = F(2);
    BW = abs(f1-f0);

    
    % get modulation function [char array]
    if iscell(P.modFn)
        modFn = P.modFn{1};
    else
        modFn = P.modFn;
    end


    % get window function for amplitude shading
    if iscell(P.winFn)
        winFn = P.winFn{1};
    else
        winFn = P.winFn;
    end
    % convert strings to function handles
    if ischar(winFn)
        winFn = str2func(winFn);
    end
    % provide fallback if invalid
    if (~isa(winFn,'function_handle') || ~exist(char(winFn),'file'))
        winFn = @rectwin;
        warning('libsproc:gen_ifpulse_multi', 'No valid windowing function found, using rectwin')
    end


    % get normalized window gain modifier
    if (length(P.gain) == 1)
        gain = P.gain;
    else
        gain = P.gain(h);
    end


    % get initial component phase
    if length(P.phase) == 1
        phase = P.phase;
    else
        phase = P.phase(h);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % define time series for current IA/IF definition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t = (0:1/fs:T)';


    % define instantaneous amplitude (IA) vector
    A = gain * winFn(length(t));


    % define instantaneous frequency (IF) vector
    switch modFn
        case 'cw'
            IF = f0*ones(size(t));
        case 'lfm'
            IF = f0*ones(size(t)) + (f1-f0)/(T).*t;
        case 'sfm'
            IF = BW./2*sin(2*pi*fMod*t) + f0+BW/2;
        case 'hfm'      % only valid for downward sweep
            IF = 1./linspace(1/f0, 1/f1, length(t))';
%        case 'spline'
%            IF = ppval(P.modFn,t);       % need to test this
        otherwise
            error('No valid modulation function found in P.modFn')
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate and add the current monocomponent pulse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % call IF approximation function
    xc = gen_ifpulse(fs, IF, phase, A);
    
    % add current pulse to time series signal
    idx = (1:length(t)) + round(t0*fs);
    x(idx) = x(idx) + xc;


end

