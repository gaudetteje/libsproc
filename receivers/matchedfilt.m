function res = matchedfilt(data,hmatch)
% MATCHEDFILT  Computes the correlation of time series data with a known signal
%
% res = matchedfilt(data,hmatch) returns the output sequence from the filter
%
% ASSUMES that data is one ROW per element (or beam)
%
% hmatch is the signal replica
%

% make the replica a column vector
hmatch = hmatch(:);

% ensure the replica is unity norm
hmatch = hmatch./norm(hmatch);

res = fftfilt(conj(flipud(hmatch)),data.').';
%res = conv(conj(flipud(hmatch)), data.').';


