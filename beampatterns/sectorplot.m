% plot sector scan of relative echo intensity, time delay, and/or sensitivity

% Basic assumptions:
%   - Transmission losses are range and frequency dependent and accounted
%     for by only spherical spreading and atmospheric absorption
%   - Transducer is a Piston-type of 4.7 cm radius, as found by Hartley and
%     Suthers (1989)
%   - Receiver(s) are omnidirectional and do not contribute to angular
%     directivity
%   - Echo strength is uniform across all frequencies

close all
clear
clc

PLOTFLAG = true;

% environmental parameters
T = 25;                 % temperature [deg. C]
rh = 50;                % relative humidity [%]
p = 101.325;            % barometric pressure (1 atm)

% user parameters
fNum = 101;
rNum = 101;
tNum = 101;
f = logspace(4,5,fNum);       % frequency points - log scale [Hz]
r = linspace(0.1,5,rNum)';         % range points - linear scale [m]
theta = linspace(-90,90,tNum);       % azimuth points - linear scale [degrees]
r1 = 1;               % TL reference range [m]

% selected region of focus
r0 = 1;                 % central region of focus in range [m]
theta0 = 0;             % central region of focus in azimuth [m]

% define region of delay sensitivity
tau = 1e-5;             % delay sensitivity for contour plot

% amplitude latency trading constant [us/dB]
alt_const = -16e-6;


%% compute range-dependent losses

% calculate range dependent spherical spreading losses (logarithmic with range)
TL_sph = 20*log10(r./r1) * ones(1,numel(f));        % TL [dB]

% calculate frequency dependent absorption and spreading losses (linear with range)
alpha = calcAbsorptionCoef(f,rh,T,p);
TL_abs = r * alpha';     % TL [dB] => [dB/m * m]

% find total combined 2-way transmission loss vs. range and frequency
TL = TL_abs + TL_sph;   % TL [dB] = dB + dB

% plot transmission loss
if PLOTFLAG
    figure;
    semilogx(r,TL)
    grid on
    %legend(num2str(1e-3*f'),'location','northwest')
    xlabel('Range (m)')
    ylabel(sprintf('Attenuation (dB re %g m)',r1))
    title('One-Way Transmission Loss in Air')
end

% convert TL to ES
TS = 0;                 % target strength [dB]
ES = -2*TL + TS;        % echo strength relative to source level [dB SPL @ r1]

% plot relative echo strength
if PLOTFLAG
    figure;
    semilogx(r,ES)
    grid on
    %legend(num2str(1e-3*f'),'location','southwest')
    xlabel('Range (m)')
    ylabel(sprintf('Relative Echo Intensity (dB re %g m)',r1))
    title(sprintf('Frequency Dependent Echo Strength for TS = %g dB',TS))
end

%% compute beam pattern directivity
D = 0.0094;             % transmit aperture
c = calcSoundSpeed(T);  % compute speed of sound in air
B = calcPistonBeam(D,f,theta(:),c);
B = 10*log10(B);

% plot beam pattern vs frequency and azimuth
if PLOTFLAG
    figure;
    plot(theta,B)
    axis([theta(1) theta(end) -30 0])
    %legend(num2str(1e-3*f'),'location','south')
    grid on
    xlabel('Azimuth (deg)')
    ylabel('Intensity (dB)')
    title(sprintf('Angular Directivity of Source with d = %g cm',D*1e2))
end


%% combine range and angular losses to form sector plot
% Z results in a 3D data cube where range, azimuth, and frequency are the
% dimensions of the NxMxL matrix
N = numel(r);
M = numel(theta);
L = numel(f);

V1 = permute(ES,[1 3 2]);
V1 = repmat(V1,[1 M 1]);

V2 = permute(B,[3 1 2]);
V2 = repmat(V2,[N 1 1]);

% combine results (by summation in dB)
V = V1 + V2;

% plot slices through volumetric data
if PLOTFLAG
    figure;
    [X,Y,Z] = meshgrid(theta,r,f*1e-3);
    xslice = [0];        % azimuth slice
    yslice = [2.5];        % range slice
    zslice = [20];      % frequency slice
    h = slice(X,Y,Z,V,xslice,yslice,zslice);
    set(gca,'clim',[-40 0])
    set(gca,'ZDir','reverse')
    set(h,'EdgeColor','none')  %'FaceColor','interp')
    xlabel('azimuth (deg)')
    ylabel('range (m)')
    zlabel('frequency (kHz)')
    title(sprintf('Relative Echo Intensity (dB re %g m)',r1))
    colorbar
%    colormap(hot)
    view(0,90)
end


%% compute the 


% find the effective transfer function for a given focus point in the range/azimuth plane
