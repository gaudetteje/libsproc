%%% Computes a ULA beam pattern by beamforming on a synthetic planar wave

clc
clear


% analysis parameters
f = 1e3;            % analysis frequency [Hz]
phi = 15;            % incident planar wave direction [degrees]
%theta = -90:0.1:90;   % angular coordinates [degrees]
theta = -180:0.1:180;   % angular coordinates [degrees] - show planar symmetry

% physical parameters
%c = 1467;           % speed of sound [m/s]
c = 1500;

% array parameters
N = 22;             % number of elements
dx = 0.15;          % element distance [m]

% array shading coefficients
%w = ones(N,1);      % rectangular window (natural array pattern)
w = chebwin(N,40);  % Dolph-Chebychev window w/ sidelobe constraints
%w = hamming(N);     % Hamming window
%w = hann(N);
w = w./sum(w);      % enforce distortionless response constraint

% plot parameters
DR = 50;            % set plot dynamic range

%%%%%%%%
%lambda = c./f;

%% construct beam


% calculate angular index
%nfft = 4096;
%L  = (-nfft/2:(nfft/2)-1);
%theta = asin(L .* c./(f.*dx.*nfft)) .* 180/pi;

% steering vectors for each look direction
delay_phi = (0:N-1)'*dx*cosd(90-theta)/c;     % delay time vector [s]
dphi = exp(1i*2*pi*f*delay_phi);       % phase steering vector to phi

% construct unit energy data vector with aperture shading weights applied
delay_theta = (0:N-1)'*dx*cosd(90-phi)/c;
dtheta = exp(1i*2*pi*f*delay_theta);
x = w .* dtheta;

%W = diag(w);        % aperture shading diagonal matrix
Y = x' * dphi;
b = db(abs(Y).^2,'power');


%% compute pattern in frequency domain
%b = abs(fftshift(fft(Y,nfft)));
%b = db(b);
%pattern = db(pattern./max(pattern));        % normalize to maximum


%

% compute half-power beamwidth (HPBW)
%HPBW = calcBeamwidth(angles,pattern)

% compute the directivity index
%DI = calcBeamDirectivity(angles,pattern)


%% plot beam pattern
figure(1)
ph = plot(theta,b);
grid on;
hold on;
xlabel('Angle (degrees)')
ylabel('Magnitude (dB)')
set(gca,'ylim',[-DR 3])
set(gca,'xlim',[theta(1) theta(end)])

color = {'b','g','r','c','m','k'};
idx = length(get(gca,'Children'));
idx = mod(idx,length(color));
set(ph,'Color',color{idx});

%%% plot beam pattern on polar coordinates
b(b < -DR) = -DR;
bpol = DR+b;

figure(2)
ph2 = polar(pi*theta/180,bpol);
hold on
set(ph2,'Color',color{idx});
view(270,90)
set(gca,'YDir','reverse')

%tilefigs(3,3)
