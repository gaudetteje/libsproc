function fc = calc_ctrfreq(x,fs)
% CALC_CTRFREQ  compute the center frequency of a time series signal
%
% f0 = calc_ctrfreq(x,fs) returns the modulation frequency (or central
%   frequency) of a signal.  f0 satisfies the intuitive criteria such that:
%       Int[0,inf] (f-f0) |X(f)|^2 df = 0
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

% default parameters
nfft = 2^12;
L = numel(x);
nfft = max(nfft,2.^nextpow2(L));

% frequency index
f = (0:nfft-1)'*fs/nfft;
f = f(1:end/2);
df = fs/nfft;

% compute spectrum
Sxx = fft(x,nfft)/sqrt(L);
Sxx = Sxx(1:end/2);

% calculate center frequncy
n = trapz(f.*(abs(Sxx).^2)) .* df;
d = trapz(abs(Sxx).^2) .* df;
fc = n./d;

%%% sanity check
res = trapz((f-fc) .* abs(Sxx).^2 .* df);
if(abs(res) > 1e5)
    warning('Residual not close to zero; fc estimate may not be accurate')
end

% compute signal energy
%E1 = trapz(abs(x).^2) ./ fs;
%E2 = trapz(abs(Sxx).^2).* df ./ fs;

if 1
    figure(1)
    plot(f*1e-3,db(abs(Sxx)))
    hold on
    grid on
    ylabel('Frequenc (kHz)')
    
end