nfft = 1024;
array = ones(1,12);

% apply aperture shading
%array = array .* hamming(length(array)).';

% random placement to reduce grating lobes

% drop channels
%array(2) = 0;%[5 12 31 40]) = 0;

% steer beam to off-boresight
%array = array .* exp(j*2*pi*1);

c = 344;
f = 20e3;%[20:10:100].* 1e3;
d = .014;

lambda = c./f;

pattern = abs(fftshift(fft(array,nfft)));
pattern = db(pattern./max(pattern));

L  = (-nfft/2:(nfft/2)-1);
angles = asin(L * c/(f*d*nfft)) * 180/pi;


figure(1)
ph = plot(angles,pattern);
grid on;
hold on;
xlabel('Angle (degrees)')
ylabel('Magnitude')


color = {'b','g','r','c','m','k'};
idx = length(get(gca,'Children'));
idx = mod(idx,length(color))+1;
set(ph,'Color',color{idx});

pattern(pattern == -Inf) = -100;
patt = 100+pattern;

figure(2)
ph2 = polar(pi*angles/180,patt);

set(ph2,'Color',color{idx});
