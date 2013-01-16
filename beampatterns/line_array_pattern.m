

% physical parameters
c = 344;            % speed of sound [m/s]

% array parameters
N = 4;              % number of elements
dx = 0.508;       % element distance [m]

% analysis parameters
theta = -90:1:90;
f = 10e3;
nfft = 1024;

%%%%%%%%
lambda = c./f;

% construct array
array = ones(1,N);
array(2) = 0;

% compute sinc pattern based on fft
pattern = abs(fftshift(fft(array,nfft)));
pattern = db(pattern./max(pattern));

% 
L  = (-nfft/2:(nfft/2)-1);
angles = asin(L .* c./(f.*d.*nfft)) .* 180/pi;


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
hold on
set(ph2,'Color',color{idx});

tilefigs(2,2)
