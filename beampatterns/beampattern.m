%close all
clc

%%%%%%%% something broke...

nfft = 2^10; %2^16;

M = 100;
array = ones(1,M);
wind = rectwin(M)';

c = 1500;
f = 50e3;
lambda = c./f;

d = lambda;
%L = (N-1)*d;       % aperture length

% generate sinc function for beam pattern
beta = abs(fftshift(fft(wind .* array,nfft)));
beta = db(beta./max(beta));

% scale
%L  = (-nfft/2:(nfft/2)-1) / nfft;
L = linspace(-1,1,nfft+1); L(end) = [];
L = L * (lambda / (d*2));

% remove beam indices outside of +/-90 degrees
idx = find(L < -1 | L > 1);
numel(idx)
L(idx) = [];
beta(idx) = [];

% replicate beam indices within +/-90 degrees


% warp angular axis by 
phi = asind(L);

% remap beam to uniform angular coordinates
phi_interp = linspace(-90,90,nfft);
beta_interp = interp1(phi,beta,phi_interp);

% beta = abs(fftshift(fft(wind .* array,nfft)));
% beta = db(beta./max(beta));
% 
% L  = (-nfft/2:(nfft/2)-1);
% phi = asin(L * c/(f*d*nfft)) * 180/pi;


%polar(phi*pi/180,B_dB)
%view(-90,90)
%set(gca,'ydir','reverse')


figure(1)
subplot(2,1,1)
ph = plot(phi,beta);
grid on;
hold on;
set(gca,'fontsize',20)
xlabel('Angle (degrees)','fontsize',20)
ylabel('Magnitude','fontsize',20)
axis([-90 90 -60 6])
set(gca,'xtick',(-90:30:90))

subplot(2,1,2)
ph = plot(phi_interp,beta_interp);
grid on;
hold on;
set(gca,'fontsize',20)
xlabel('Angle (degrees)','fontsize',20)
ylabel('Magnitude','fontsize',20)
axis([-90 90 -60 6])
set(gca,'xtick',(-90:30:90))

%text(15,-5,'\beta = 1.5\circ','fontsize',22,'fontweight','bold')
%hold on;
%plot(phi_interp, beta_interp, '--r')


%color = {'b','g','r','c','m','k'};
%idx = length(get(gca,'Children'));
%idx = mod(idx,length(color))+1;
%set(ph,'Color',color{idx});
% 
% pattern(pattern == -Inf) = -100;
% patt = 100+pattern;
% 
% figure(2)
% ph2 = polar(pi*angles/180,patt);
% view(-90,90)

%set(ph2,'Color',color{idx});

