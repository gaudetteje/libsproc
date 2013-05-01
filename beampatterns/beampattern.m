close all

nfft = 2^16;

% random placement to reduce grating lobes

% drop channels
%array(2) = 0;%[5 12 31 40]) = 0;

% steer beam to off-boresight
%array = array .* exp(j*2*pi*1);

c = 1500;
f = [20:10:100].* 1e3;
lambda = c./f;

d = .02;%lambda/2;%.014;

M = 22;

array = ones(1,M);
wind = rectwin(M)';

pattern = abs(fftshift(fft(wind .* array,nfft)));
pattern = db(pattern./max(pattern));

L  = (-nfft/2:(nfft/2)-1);
angles = asin(L * c/(f*d*nfft)) * 180/pi;


figure(1)
ph = plot(angles,pattern);
grid on;
hold on;
set(gca,'fontsize',20)
xlabel('Angle (degrees)','fontsize',20)
ylabel('Magnitude','fontsize',20)
axis([-90 90 -60 6])
set(gca,'xtick',[-90:30:90])

%text(15,-5,'\beta = 1.5\circ','fontsize',22,'fontweight','bold')


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
