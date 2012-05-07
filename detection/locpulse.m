function pulses=locpulse(s, fs, fc, varargin)
% LOCPULSE  This program detects the sample numbers and times for CW sonar
% pulses.
%
% PULSES = LOCPULSE(T,S,FC,THRES) uses a matched filter for a sinusoid of
% frequency FC to detect pulses in the signal vector S.  PULSES is a
% structure of starting values where PULSES.index contains the sample
% indices and PULSES.time contains the times based on the time vector T.

% Author:   Jason Gaudette
% Company:  Naval Undersea Warfare Center (Newport, RI)
% Phone:    401.832.6601
% Email:    gaudetteje@npt.nuwc.navy.mil
% Date:     20061002
%

hlen = 0.001;
if nargin > 3
    hlen = varargin{1};
end

% generate matched signal
t_h = [0 : 1/fs : hlen];
h = sin(2*pi*fc*t_h);

tmin = length(h);
if nargin > 4
    tmin = varargin{2};
end

% run matched filter over data
y = matchedfilt(s,h);

gamma = 0.95 * max(y);  % assumes all pings are same amplitude per plot
res = find(y > gamma) - length(t_h);

smin = tmin * fs;
pulses = res(find(diff(res) > smin)+1);
