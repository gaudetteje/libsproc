function [DI, D] = calcBeamDirectivity(phi,pattern)
% CALCBEAMDIRECTIVITY  computes the directivity index of a 1D beam pattern
%
%

D = inf;
DI = 10*log10(D);