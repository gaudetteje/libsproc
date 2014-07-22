function beta = calc_rmsband(x,fs)
% CALC_RMSBAND  compute the center frequency of a time series signal
% (based on Rhiaczek)

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

% calc center frequency
fc = calc_ctrfreq(x,fs);

% calculate rms bandwidth
n = trapz((f-fc).^2 .* (abs(Sxx).^2)) .* df;
d = trapz(abs(Sxx).^2) .* df;
beta = 2*pi*sqrt(n./d);

plot(f,abs(Sxx))
grid on
