
% physical parameters
c = 1467;           % speed of sound [m/s]
%c = 1500;

% array parameters
N = 11;             % number of elements
dx = 0.15;          % element distance [m]

% array shading coefficients
%w = ones(N,1);      % rectangular window (natural array pattern)
w = chebwin(N,30);  % Dolph-Chebychev window w/ sidelobe constraints
w = w./sum(w);      % enforce distortionless response constraint

% analysis parameters
theta = -90:1:90;
f = 5e3;
nfft = 4096;

% plot parameters
DR = 40;            % set plot dynamic range
phi = 60;            % "look" or steering angle

%%%%%%%%
lambda = c./f;
W = diag(w);

% construct beam
%d = ones(1,N);      % steering vector to broadside (natural response)
delta = (0:N-1)*dx*cosd(90-phi)/c;     % delay vector
s = exp(1i*2*pi*f*delta);       % steering vector to phi

Y = s * W;%* ones(N,1);
%Y = ones(1,N) .* w';



%% compute pattern in frequency domain
pattern = abs(fftshift(fft(Y,nfft)));
pattern = db(pattern);
%pattern = db(pattern./max(pattern));        % normalize to maximum

% calculate angular index
L  = (-nfft/2:(nfft/2)-1);
angles = asin(L .* c./(f.*dx.*nfft)) .* 180/pi;

%

% compute half-power beamwidth (HPBW)
%HPBW = calcBeamwidth(angles,pattern)

% compute the directivity index
%DI = calcBeamDirectivity(angles,pattern)


%%% plot beam pattern on 
figure(1)
ph = plot(angles,pattern);
grid on;
hold on;
xlabel('Angle (degrees)')
ylabel('Magnitude (dB)')
set(gca,'ylim',[-DR 3])

color = {'b','g','r','c','m','k'};
idx = length(get(gca,'Children'));
idx = mod(idx,length(color));
set(ph,'Color',color{idx});

%%% plot beam pattern on polar coordinates
pattern(pattern < -DR) = -DR;
patt = DR+pattern;

figure(2)
ph2 = polar(pi*angles/180,patt);
hold on
set(ph2,'Color',color{idx});
view(270,90)
set(gca,'YDir','reverse')

%tilefigs(3,3)
