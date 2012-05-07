function w = raisedcos(n,varargin)
% RAISEDCOS  returns a raised cosine window function given number of points, n,
% and percent taper, p.
%
% w = raisedcos(100,3) returns a 100 point column vector, w, with 3% taper
% applied to the front and back
%
% w = raisedcos(100) assumes a 6.25% taper
%

p = 6.25;
if (nargin == 2); p = varargin{1}; end;

w = ones(n,1);

% find the p-percent index into n-point array
idx = round(n*p/100);
if ~idx; error('Percent taper too large or length too small'); end

% apply front taper
w(1:idx+1) = w(1:idx+1) .* cos(pi/2*(0:1/idx:1) - pi/2)';

% apply back taper
w(end-idx:end) = w(end-idx:end) .* cos(pi/2*(0:1/idx:1))';
