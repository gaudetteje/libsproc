function y = rms(x)
% RMS  calculates the root mean square value of an input sequence

y = sqrt(mean(x.^2));