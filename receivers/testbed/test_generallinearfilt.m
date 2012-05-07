%% calculate and apply GLF/VRDR to signal
close all

%% FM1
fs = 250e3;%1e6;
T = 0.002;
f0 = 50e3;
f1 = 25e3;
phi0 = 0;

t = (0:1/fs:T)';
A = .5*raisedcos(length(t),20);

B = abs(f1-f0);
a = T*(f0*f1)/B;
b = T*f1/B;
IF1 = a./(t+b);

hmatch1 = real(gen_ifpulse(fs,IF1,phi0,A));
hmatch1 = hmatch1 ./ max(hmatch1);

%% FM2
f0 = 100e3;
f1 = 50e3;

B = abs(f1-f0);
a = T*(f0*f1)/B;
b = T*f1/B;
IF2 = a./(t+b);

hmatch2 = real(gen_ifpulse(fs,IF2,phi0,A));
hmatch2 = hmatch2 ./ max(hmatch2);

%% FM1+2
hmatch = (hmatch1+hmatch2);
hmatch = hmatch ./ max(hmatch);


%% pulse train
idx = 5;
L = 6e3;
time = (0:L-1)'./fs;
data = zeros(L,1);
data(1000:3000) = hmatch;
data(3000:5000) = hmatch;
data((3000:5000)+idx) = data((3000:5000)+idx) + hmatch;

data = data + .1*randn(size(data));


%% matched filter data
figure(1);
subplot(2,1,1)
plot(time*1e3,data)
axis([0 10 -2 2])
grid
Y = matchedfilt(data,hmatch);
subplot(2,1,2)
plot((time-.002)*1e3,Y)
axis([0 10 -20 20])
grid


% matched filter
alpha = 1;
YY = generallinearfilt(data,hmatch,alpha);

time = time-1/fs;
hold on
plot((time-.002)*1e3,YY,'r')
hold off


%% autocorr
figure(2);
[xc1,lag] = xcorr(hmatch1,hmatch1);
xc2 = xcorr(hmatch2,hmatch2);
xc = xcorr(hmatch,hmatch);

xc1 = xc1./max(xc1);
xc2 = xc2./max(xc2);
xc = xc./max(xc);

tl = 1e6*lag./fs;
xc1env = abs(hilbert(xc1));
xc2env = abs(hilbert(xc2));
xcenv = abs(hilbert(xc));

subplot(3,1,1)
plot(tl,xc1,'b','linewidth',2)
hold on; grid on;
plot(tl,xc1env,'b')
title('Autocorrelation - FM1')
axis([-200 200 -1.2 1.2])

subplot(3,1,2)
plot(tl,xc2,'r','linewidth',2)
hold on; grid on;
plot(tl,xc2env,'r')
title('Autocorrelation - FM2')
axis([-200 200 -1.2 1.2])

subplot(3,1,3)
plot(tl,xc,'k','linewidth',2)
hold on; grid on;
plot(tl,xcenv,'k')
title('Autocorrelation - FM1+2')
axis([-200 200 -1.2 1.2])
%legend('FM1','FM1 envelope','FM2','FM2 envelope')  %,'FM1+2')
xlabel('Time (\mus)')