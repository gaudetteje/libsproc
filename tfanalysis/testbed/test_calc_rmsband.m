clc
clear
close all

% define simulated parameters
T = .01;%.01;
pLen = .001; %.005;
f0 = 100e3 + 10e3*randn;     % center frequency
fsrange = logspace(log10(200e3),log10(1e6),20);

% define test signal
P.time = [0 pLen]';
P.freq = [f0; f0];
P.modFn = 'cw';
P.winFn = 'hamming'; %'raisedcos';
P.gain = 1;
P.phase = 0;


%% iterate over sampling rate to confirm no bias error exists

% init arrays
rms_est = [];
fc_est = [];
beta_est = [];

for fs = fsrange;
    % generate time series vector
    t = (0:1/fs:T)';
    
    % construct a random walk (to add bandwidth)
%     N = 1:numel(t);
%     f_i = fc + cumsum(randn(N,1));

    % generate test signal
    x = gen_ifpulse_multi(fs, P);
    L = (T-pLen).*fs;
    x = [zeros(round(L/2),1); real(x); zeros(round(L/2),1)];    % zero pad

    % confirm center frequency and rms bandwidth calculations
    [beta,fc] = calc_rmsband(x,fs);

%    rms_est(end+1) = rms(x);
    fc_est(end+1) = fc;
    beta_est(end+1) = beta;
end


figure
subplot(2,1,1)
plot(t*1e3,x)
xlabel('Time (ms)')
ylabel('Amplitude')
title(sprintf('(f_c = %g kHz, pLen = %g s)',fc*1e-3,pLen))

subplot(2,1,2)
nfft = 2^(nextpow2(numel(x)));
spectrogram(x,L,round(L*.9),nfft,fs,'yaxis')
clim = get(gca,'clim');
set(gca,'clim',[clim(2)-40 clim(2)]);
title(sprintf('NFFT = %d',nfft))


figure
subplot(3,1,1)
plot(fsrange*1e-3,fc_est*1e-3,'.')
hold on
plot([fsrange(end)]*1e-3,[f0 f0]*1e-3,'--r')
ylabel('F_c (kHz)')
title(sprintf('(f_c = %g kHz, pLen = %g ms)',fc*1e-3,pLen*1e3))

subplot(3,1,2)
plot(fsrange*1e-3,beta_est*1e-3,'.')
hold on
%plot([0 fsrange(end)]*1e-3,1e-3*[beta_hat beta_hat],'--r')
ylabel('Beta (kHz)')

subplot(3,1,3)
plot(fsrange*1e-3,fc_est./beta_est,'.')
hold on
xlabel('Sampling Rate (kHz)')
ylabel('Q (F_c / Beta)')
