function res = inversefilt(data,hmatch)
% INVERSEFILT  Computes the correlation of time series data with a known signal
%
% res = inversefilt(data,hmatch) returns the output sequence from the filter
%
% ASSUMES that data is one ROW per element (or beam)
%
% hmatch is the signal replica
%

% make the replica a column vector
hmatch = hmatch(:);

% ensure the replica is unity norm
hmatch = hmatch./norm(hmatch);

res = filter(1, flipud(hmatch), data.').';         % inverse
