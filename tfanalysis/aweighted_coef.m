function [b,a] = aweighted_coef(fs)
% FILT_AWEIGHTED designs filter coefficients for an a-weighted filter
%
% [b,a] = aweighted_coef(fs) designs an analog filter and transforms to a
% digital filter using the bilinear transform

%Analog A-weighting filter according to IEC/CD 1672.
f1 = 20.598997; 
f2 = 107.65265;
f3 = 737.86223;
f4 = 12194.217;
A1000 = 1.9997;
pi = 3.14159265358979;
NUM = [ (2*pi*f4)^2*(10^(A1000/20)) 0 0 0 0 ];
DEN = conv([1 +4*pi*f4 (2*pi*f4)^2],[1 +4*pi*f1 (2*pi*f1)^2]); 
DEN = conv(conv(DEN,[1 2*pi*f3]),[1 2*pi*f2]);

%Bilinear transformation of analog design to get the digital filter. 
[b,a] = bilinear(NUM,DEN,fs);
