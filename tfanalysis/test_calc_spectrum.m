clear
%clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% signal parameters
ts.fs = 2^18;%250e3;%               % Sampling frequency
L = ts.fs;                          % Length of signal
N = 2;                              % Number of channels
A = 10*sqrt(2);
f0 = 2^13;%25e3;%
ts.time = (0:L-1).'./ts.fs;         % time vector
x = A.*cos(2*pi*f0*ts.time);        % create pure tone



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create band limited AWGN
randn('state',0);
w = 1e-7 * randn(L,N);
% [N,Wn] = buttord([60e3 80e3]*2/ts.fs,[50e3 100e3]*2/ts.fs,0.1,90);
% [b,a] = butter(N,Wn);
% w = filter(b,a,w);
ts.data = x * ones(1,N) + w;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fname = '0dBV_awg.aam';%'noisefloor.aam';
% ts = ap_read_wave(fname);
% ts.fs = ts.sample_rate;
fprintf('\n\n\nVrms = %g (%gdB)\n',rms(ts.data),db(rms(ts.data)))

% spectral parameters
nfft = 2^14;%32768;
win = 'hanning';%'hamming';%'blackmanharris';%'flattopwin';%'rectwin';%
%%%%'blackman7';%
band = 'half';%'full';%
navg = 4;
overlap = 0.5;
acdc = 'ac';
units = 'Vrms^2';%Vrms/rthz';

% calculate spectrum
fd = calc_spectrum(ts,nfft,win,navg,overlap,band,acdc,units);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot spectrum
figure(10)
plotnum = length(get(gca,'Children'));
colorlist = get(gca,'ColorOrder');
numcolors = size(colorlist,1);
nextcolor = colorlist(mod(plotnum, numcolors)+1,:);

% plot results
plot(fd.freq,fd.magdb,'Color',nextcolor)
title(sprintf('Fs=%g - %d point FFT (%s, %s)',ts.fs,nfft,win,band))
xlabel('Frequency (Hz)')
ylabel(['dB' fd.units])
%ylabel(sprintf('Magnitude (dBVpk/rt(%g)',fd.resolution))
grid on
%axis([-fd.fs/2 fd.fs/2 -60 10])
%axis([0 120000 -100 -20])
%axis([4000 6000 -140 20])
hold on
legend show

frange = [60e3 80e3];
totnoise = total_spectral_noise(fd,frange);
fprintf('Total RMS noise in band %g to %g is %gV (%g dBV)\n',frange(1),frange(2),totnoise,db(totnoise))

fd
fd.win
