function [x,phi] = gen_ifpulse(fs,IF,varargin)
% GEN_IFPULSE  generate modulated signals given desired instantaneous frequency
%
% x = gen_ifpulse(fs,IF) defines a time series waveform, x, based on the
%       instantaneous frequency definition, IF, and sampling rate, fs
%
% x = gen_ifpulse(fs,IF,phi0,IA) also specifies optional parameters for
%       initial phase constant, phi0, and instantantaneous amplitude, IA.
%       IA may be either a scalar constant or a time-dependent vector.
%
% [x,phi] = gen_ifpulse(...) also returns the phase function, phi
%
% Inputs:
%   fs   -  sampling rate [Hz]
%   IF   -  desired instantaneous frequency vector [Hz]
%   phi0 -  initial phase [radians]
%   a    -  amplitude vector (scalar or length(IF))
%
% Outputs:
%   x    -  analytic (complex) phase modulated signal
%   phi  -  calculated phi(t) based on numerical integration of IF
%
% This function generates monocomponent pulses based on the desired instantaneous
% frequency (IF).  Any valid function of frequency vs. time, IF(t), can be
% converted to a frequency modulated pulse by applying phase modulation.
% Multicomponent signals can be constructed by combining separate monocomponent
% pulses.
%
% Method:
%   x(t) = A(t) e^(j phi(t))
%
% where:
%   phi(t) = 2 PI Int{ IF(t) dt }
%
%
% Usage Examples:
%
% Define basic parameters
%    fs = 1e6;                      % sampling rate [Hz]
%    T = 1e-3;                      % pulse length [seconds]
%    t = (0:1/fs:T).';              % time reference vector
%    L = numel(t);                  % pulse length [samples]
%
% Example 1 - CW pulse
%    fc = 50e3;                     % pulse center frequency
%    IF = fc*ones(L,1);
%
% Example 2 - LFM
%    f0 = 50e3;                     % start frequency
%    f1 = 100e3;                    % stop frequency
%    IF = f0*ones(L,1) + (f1-f0)/(T).*t;
%
% Example 3 - SFM
%    fc = 50e3;                     % center frequency
%    B = 60e3;                      % bandwidth
%    fm = 5000;                     % sinusoidal modulation rate
%    IF = B./2*sin(2*pi*fm*t) + fc;
%
% Example 4 - HFM
%    f0 = 100e3;                    % start frequency
%    f1 = 50e3;                     % stop frequency
%    B = abs(f1-f0);                % bandwidth
%    a = T*(f0*f1)/B;
%    b = T*f1/B;
%    IF = a./(t+b);
%
% Generate the waveform
%    IA = raisedcos(L);             % instantaneous amplitude
%    phi0 = pi/2;                   % initial phase constant
%    x = gen_ifpulse(fs,IF,phi0,IA);
%
% see also gen_ifpulse_multi, gen_pulsetrain

% check for errors in definition
if any(isinf(IF)), error('Cannot numerically integrate instantaneous frequency.  Vector contains Inf values.'), end

% set default parameters
nArgs = length(varargin);
optargs = {-pi/2 1};

% assign optional input parameters
optargs(1:nArgs) = varargin;
[phi0, a] = optargs{:};

% integrate the function, IF(t), to find phi(t)
phi = cumtrapz(IF)./fs;

% generate FM waveform
x = a .* exp(1i*2*pi*phi + 1i*phi0);

