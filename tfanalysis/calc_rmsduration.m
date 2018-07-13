function tau = calc_rmsduration(x,fs)
% CALC_RMSDURATION  computes the RMS signal duration
%
% Reference:
%     Rihaczek, A. (1996). Principles of High-Resolution Radar. Norwood, 
%     MA: Artech House.

% normalize sampling rate if none entered
if nargin == 1
    fs = 1;
end

% force data to a column vector
x = x(:);

dt = 1/fs;
t = (1:numel(x))'.*dt;

% calc RMS duration
n = trapz(t.^2 .* abs(x).^2)*dt;
d = trapz(abs(x).^2)*dt;

tau = n/d;
