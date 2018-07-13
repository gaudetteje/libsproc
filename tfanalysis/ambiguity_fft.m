function [Fu_up, Fu_down, Fv_up, Fv_down, X_up, X_down] = ambiguity_fft(fc, bands, replicasAB, LFM_data, header)
% AMBIGUITY_FFT calculates the Ambiguity Function via FFT method using desired Doppler shift fc
%
% Must load saved ERAS processing data into Workspace 
%  For example:  ...\C821 Work\ERAS_SAS\Dodge_Pond_Data_Processing\LFM_HFM_Processing_Code\LFM_1_1.mat
%
% Ambiguity function is calculated using individual range cuts at various Doppler shift values.
% We let u(t) be the baseband transmited signal and v(t) be the signal with
% which the matched filter is matched(replica).  The algorithm to generate
% each range cut, |X(Tau,f)|, where Tau is the time delay and f is the
% Doppler shift frequency:
% (a) Find FFT{u(t)exp(j*2*pi*f*t)} = Fu(theta) ----> FFT of Doppler-shifted baseband signal
% (b) Find FFT{v(t)} = Fv(theta) ----> FFT of replica
% (c) Find Fx(theta) = Fu(theta)Fv*(theta) = FFT{X(Tau,f)} ----> Product of FFT of signal with the complex conjugate of FFT of replica
% (d) Find X(Tau,f) = IFFT{Fx(theta)}
% (e) Take the absolute value of X(Tau,f) to find Ambiguity Function, |X(Tau,f)|
%
% Define a time array, t, of length N.  
% N must be an integer power of 2 (N=2^M where M is a positive integer) and extends from 0 to T.
% T must be equal to or great than the sum of durations of u(t) and v(t)
% N must be chosen such that N/T >2F where F is the highest frequency of u(t) and v(t)
% The larger the N, the better the range cut will look
%
% Code based on radar systems course notes:
%   http://www.ece.uah.edu/courses/material/EE619/index.htm

T = 2*(length(bands(1).baseband_data) * 1/(header.fs/47));
%T = 2*(length(bands(1).baseband_data) * 1/header.fs);
%N = 2048;
N = length(bands(1).baseband_data); % 433
t = (0:T/N:T-T/N); % [1x433]

%NFFT = N;
NFFT=2048;

% Apply Doppler Shift to Up and Down data
% size(bands(1).baseband_data) = [377x112]
Doppler_shift = exp(1i*2*pi*fc*t).'; % [433x1]

% Zero-pad

% ZP_up = zeros((N-(length(bands(1).baseband_data))), 112);
% ZP_down = zeros((N-(length(bands(2).baseband_data))),112);
% data_up = cat(1,bands(1).baseband_data, ZP_up);
% data_down = cat(1,bands(2).baseband_data, ZP_down);
data_up = bands(1).baseband_data; %[433x112]
data_down = bands(2).baseband_data;

DS_up = zeros(size(data_up)); % [433x112]
DS_down = zeros(size(data_down));
for i=1:112
    DS_up(:,i) = (Doppler_shift).*(data_up(:,i));
    DS_down(:,i) = (Doppler_shift).*(data_down(:,i));
end



% Compute FFT of Doppler-shifted data
Fu_up = fftshift(fft(sum((DS_up(:,:)).'),NFFT)); % [1x2048]
Fu_down = fftshift(fft(sum((DS_down(:,:)).'),NFFT));

% % Plot FFTs of Doppler-shifted data 
% figure; plot((-NFFT/2:NFFT/2-1)/NFFT*header.fs, 20*log10(abs(Fu_up)))
% grid on;
% 
% figure; plot((-NFFT/2:NFFT/2-1)/NFFT*header.fs, 20*log10(abs(Fu_down)))
% grid on;


% Compute FFT of replicas
% Fv_up = 20*log10(abs(fftshift(fft(sum((replicasAB(1).baseband_data(:,:)).'),NFFT))));
% Fv_down = 20*log10(abs(fftshift(fft(sum((replicasAB(2).baseband_data(:,:)).'),NFFT)))); 
Fv_up = fftshift(fft(((replicasAB(1).baseband_data(:,:)).'),NFFT)); % [1x2048]
Fv_down = fftshift(fft(((replicasAB(2).baseband_data(:,:)).'),NFFT)); 

% % Plot FFTs of replicas 
% figure; plot((-NFFT/2:NFFT/2-1)/NFFT*header.fs, 20*log10(abs(Fv_up)))
% grid on;
% 
% figure; plot((-NFFT/2:NFFT/2-1)/NFFT*header.fs, 20*log10(abs(Fv_down)))
% grid on;

% Compute Ambiguity Function  |Chi(Tau',fc)| = X = |IFFT(Fu*Fv')| 
X_up = abs(fftshift(ifft(Fu_up.*conj(Fv_up)))); % [1x2048]
X_down = abs(fftshift(ifft(Fu_down.*conj(Fv_down))));

% Plot Range cut of the Ambiguity Function at a Doppler shift of fc
% figure; 
% %hold on;
% plot(X_up, 'k')
% hold on;
% plot(X_down, 'r')
% 
% figure; 
% %hold on;
% plot(-T/2:T/N:(T/2-T/N), X_up, 'k')
% hold on;
% plot(-T/2:T/N:(T/2-T/N), X_down, 'r')
% 
% figure; 
% %hold on;
% plot(1415*(-T/2:T/N:(T/2-T/N))/2, X_up, 'k')
% hold on;
% plot(1415*(-T/2:T/N:(T/2-T/N))/2, X_down, 'r')
% 
% figure; 
% %hold on;
% plot(1415*t/2, X_up, 'k')
% hold on;
% plot(1415*t/2, X_down, 'r')
% 
% figure; 
% %hold on;
% plot(1415*t/2-53.6285, X_up, 'k')
% hold on;
% plot(1415*t/2-53.6285, X_down, 'r')

% tt = [-T/2:T/N:T/2-T/N];
% %figure; 
% %hold on;
% plot(1415*tt/2, X_up, 'blue')
% hold on;
% plot(1415*tt/2, X_down, 'magenta')
