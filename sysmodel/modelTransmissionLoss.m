close all
clear
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% define user parameters and desired system response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data params
fs = 192e3;             % sampling rate [Hz] - should match your data!
fInt = 1e3;             % frequency domain sampling interval [Hz]

% environmental params
h_r = 20;               % relative humidity [%]
T = 20;                 % ambient temperature [deg C]
P = 101.325;            % barometric pressure [kPa]
range = 2;              % range of acoustic transmission [m]

% system model params
Nb = 10;                % number of zeros
Na = 10;                % number of poles
maxIter = 10;           % maximum iterations with Steiglitz-McBride algorithm

% plot params
cLim = [-80 -40];


% calculate transmission loss due to absorption at specified range
f = (0:fInt:fs/2);
calcAbsorptionCoef(f, h_r, T, P);       % plot first
alpha = calcAbsorptionCoef(f, h_r, T, P);
TL = 20*log10(range) + alpha.*range;     % assume spherical (free-field) propagation


% define desired magnitude response from calculated transmission loss
Hd = 10.^(-TL/20);
Hd = [Hd; flipud(Hd(2:end-1))];      % mirror full spectrum (removing point at Nyquist)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% estimate system transfer function poles/zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Hd = abs(Hd(:));               % ensure magnitude only; convert to column vector
L = length(Hd);
W = dftmtx(L); Wb = W; Wa = W;
Wb(:,Nb+2:L) = []; Wa(:,Na+2:L) = [];

% generate the autocorrelation function
r = ifft(Hd.^2);

% construct an initial system model by Levinson-Durbin (AR), follow with Prony (ARMA)
aL = levinson(r,floor(L/2));
hL = impz(1,aL,Nb+2*Na+2);
[b,a] = prony(hL,Nb,Na);

% iteratively refine pole/zero positions with frequency domain Steiglitz-McBride (ARMA)
for i = 1:maxIter,
    [Hi,w] = freqz(b,a,L,'whole');
    Hai = freqz(1,a,L,'whole');
    Pi = exp(1i*angle(Hi));
    HdPi = Hd.*Pi;
    b = (diag(Hai)*Wb)\HdPi; B = fft(b,L);
    a = (diag(HdPi.*Hai)*Wa)\(diag(Hai)*Wb*b);
end

% force real filter coefficients (Hd should be symmetric)
if (sum(imag(b)) + sum(imag(a)) > 1e-10)
    warning('Poles and/or zeros not entirely real.  Possibly throwing away significant imaginary parts.')
    fprintf('Note:  The desired magnitude response, Hd, should be kept symmetric to ensure real coefficients!\n')
end
b = real(b.');
a = real(a.');

% scale all coefficients evenly, forcing a0=1
b = b/a(1);
a = a/a(1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% post-synthesis model analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
[H,w4] = freqz(b,a,4*L,'whole');        % interpolated frequency response
plot(w4/pi,abs(H),'k-',w/pi,Hd,'k--')
grid on;
title('ARMA System Model, H(z) = S(z) = B(z)/A(z)')
xlabel('Normalized frequency (w/pi)')
ylabel('Magnitude')
legend('Model magnitude','Desired magnitude','location','Best')

H = H(1:4:4*L);     % downsample frequency domain

figure
plot(w/pi,Hd - abs(H),'r')
grid on;
title('ARMA System Model Error, Hd(z) - H(d)')
xlabel('Normalized frequency (w/pi)')
ylabel('Magnitude')

figure;
zplane(b,a);
title('Pole/zero plot')

fprintf('Sum of sq. error:  %g\n', sum((Hd - abs(H)).^2))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test model using a synthetic wideband signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = (1:1024)./fs;
s = zeros(1024,1);
s_hat = fmlin(512);
s(257:768) = real(s_hat);
n = 0.001 .* randn(1024,1);
x = s + n;

y = filter(b,a,x);

figure;
plot(t,x,'b',t,y,'g')
grid on;
ylabel('Amplitude')
xlabel('Time (sec)')
title('Time series waveforms')
legend('Input signal','Output signal')

figure;
pwelch(x)
hold on;
pwelch(y)
grid on;
legend('Input signal','Output signal')

figure;
spectrogram(x,64,60,[],fs,'yaxis')
caxis(cLim)

figure;
spectrogram(y,64,60,[],fs,'yaxis')
caxis(cLim)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reverse propagation losses using the inverse filter, H(z) = 1/S(z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
H_inv = freqz(a,b,4*L,'whole');
plot(w4/pi,abs(H_inv),'k-',w/pi,1./Hd,'k--')
grid on;
title('ARMA Inverse System Model, H(z) = 1/S(z) = A(z)/B(z)')
xlabel('Normalized frequency (w/pi)')
ylabel('Magnitude')
legend('Model magnitude','Desired magnitude','location','Best')

H_inv = H_inv(1:4:4*L);

figure
plot(w/pi,1./Hd - abs(H_inv),'r')
grid on;
title('ARMA Inverse System Model Error, Hd(z) - H(z)')
xlabel('Normalized frequency (w/pi)')
ylabel('Magnitude')

x_hat = filter(a,b,y);

figure;
plot(t,x,'b',t,x_hat,'r--')
grid on;
ylabel('Amplitude')
xlabel('Time (sec)')
title('Time series waveforms')
legend('Original signal','Reconstructed signal')

figure;
plot(t,x - x_hat,'-k')
ylabel('Amplitude')
xlabel('Time (sec)')
grid on;
title('Waveform reconstruction error')
