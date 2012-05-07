clear
clc
close all

fs = 192e3;

% define desired magnitude response
f = (0:1e3:fs);
Hd = [ones(1,length(f)/4-25) (25:-1:1)./25 (0:25)./25 ones(1,length(f)/4-25)];
%Hd = 1:length(f)/2+1;
%Hd = 10.*ones(1,length(f)/2+1);
Hd = [Hd fliplr(Hd(2:end))];

plot(f,Hd)
grid on


[b,a] = iirmag(Hd,20,20);
b = real(b);
a = real(a);


% post analysis
figure;
[H,W] = freqz(b,a,f,fs);
plot(abs(H))
grid

figure;
plot(Hd - abs(H))
grid
