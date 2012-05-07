function w = blackman7(n)
% Compute the 7-Term Blackman-Harris Window for N-Points
%
%	        Input Arguments:
%
%			n   => Window Length;
%
%      		Output Arguments:
%
%			w   => Window Coefficients;

%		References:								
%                               [1] O. M. Solomon, Jr., "The use of DFT windows in 
%                                   signal-to-noise ratio and harmonic distortion 
%                                   computations", IEEE Transactions on Instrumentation 
%                                   and Measurement, vol. 43, no. 2, pp 194-199, April 
%                                   1994.
%
%       			[2] F. J. Harris, "On the use of windows for harmonic	
%				    analysis with the Discrete Fourier Transform",	
%				    Proceedings of the IEEE, Vol. 66, No. 1, Jan 1978.	
%											
%				[3] A. H. Nuttal, "Some windows with very good sidelobe	
%				    behavior", IEEE Trans. on Acoustics, Speech & Signal
%				    Processing, Vol. ASSP-29, No. 1, Feb 1981.		



a(1) = 0.271051400693424;
a(2) = 0.433297939234485;
a(3) = 0.218122999543110;
a(4) = 0.065925446388031;
a(5) = 0.010811742098371;
a(6) = 0.000776584825226;
a(7) = 0.000013887217352;

%w = ((a(1) - a(2)*cos(2*pi*(0:n-1)/n)) + a(3)*cos(4*pi*(0:n-1)/n)) - a(4)*cos(6*pi*(0:n-1)/n);
%w = ((w + a(5)*cos(8*pi*(0:n-1)/n)) - a(6)*cos(10*pi*(0:n-1)/n)) + a(7)*cos(12*pi*(0:n-1)/n);

w = (a(7)*cos(12*pi*(0:n-1)/n) - a(6)*cos(10*pi*(0:n-1)/n)) + a(5)*cos(8*pi*(0:n-1)/n);
w = ((w - a(4)*cos(6*pi*(0:n-1)/n)) + a(3)*cos(4*pi*(0:n-1)/n)) - a(2)*cos(2*pi*(0:n-1)/n);
w = w' + a(1);
