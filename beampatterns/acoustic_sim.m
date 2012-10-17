function acoustic_sim
% Plot sector scan of relative echo intensity, time delay, and/or sensitivity

% Basic assumptions:
%   - Transmission losses are range and frequency dependent and accounted
%     for by only spherical spreading and atmospheric absorption
%   - Transducer is a Piston-type of 4.7 cm radius, as found by Hartley and
%     Suthers (1989)
%   - Receiver(s) are omnidirectional and do not contribute to angular
%     directivity
%   - Transmitted and reflected echo strength is uniform across all frequencies

close all
clear
clc

PLOTFLAG = true;

set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)


% physical parameters
D = 0.0094;             % transmit piston aperture

% environmental parameters
T = 25;                 % temperature [deg. C]
rh = 50;                % relative humidity [%]
p = 101.325;            % barometric pressure (1 atm)

% plotting parameters
fNum = 101;
rNum = 101;
tNum = 101;
f = logspace(4,5,fNum);         % frequency points - log scale [Hz]
rho = linspace(0.1,5,rNum)';   % range points - linear scale [m]
theta = linspace(-90,90,tNum);  % azimuth points - linear scale [degrees]
rhoRef = .1;                   % TL reference range [m]

% derive parameters for later use
N = numel(rho);
M = numel(theta);
L = numel(f);


%% compute range-dependent losses

% calculate range dependent spherical spreading losses (logarithmic with range)
TL_sph = 20*log10(rho./rhoRef) * ones(1,numel(f));        % TL [dB]

% calculate frequency dependent absorption and spreading losses (linear with range)
alpha = calcAbsorptionCoef(f,rh,T,p);
TL_abs = rho * alpha';     % TL [dB] => [dB/m * m]

% find total combined 2-way transmission loss vs. range and frequency
TL = TL_abs + TL_sph;   % TL [dB] = dB + dB

% plot transmission loss
if PLOTFLAG
    figure;
    semilogx(rho,TL)
    grid on
    xlabel('Range (m)')
    ylabel(sprintf('Attenuation (dB re %g m)',rhoRef))
    title('One-Way Transmission Loss in Air')
    if numel(f) <= 11
        legend([num2str(round(1e-3*f')) ones(L,1)*' kHz'],'location','northwest')
    end
end

% convert TL to ES
TS = 0;                 % target strength [dB]
ES = -2*TL + TS;        % echo strength relative to source level [dB SPL @ r1]

% plot relative echo strength
if PLOTFLAG
    figure;
    semilogx(rho,ES)
    grid on
    xlabel('Range (m)')
    ylabel(sprintf('Relative Echo Intensity (dB re %g m)',rhoRef))
    title(sprintf('Frequency Dependent Echo Strength for TS = %g dB',TS))
    if numel(f) <= 11
        legend([num2str(round(1e-3*f')) ones(L,1)*' kHz'],'location','southwest')
    end
end

%% compute beam pattern directivity
c = calcSoundSpeed(T);  % compute speed of sound in air
B = calcPistonBeam(D,f,theta(:),c);
B = 10*log10(B);

% plot beam pattern vs frequency and azimuth
if PLOTFLAG
    figure;
    plot(theta,B)
    axis([theta(1) theta(end) -30 0])
    grid on
    xlabel('Azimuth (deg)')
    ylabel('Intensity (dB)')
    title(sprintf('Angular Directivity of Source with d = %g cm',D*1e2))
    if numel(f) <= 11
        legend([num2str(round(1e-3*f')) ones(L,1)*' kHz'],'location','south')
    end
end


%% combine range and angular losses to form sector plot
% Z results in a 3D data cube where range, azimuth, and frequency are the
% dimensions of the NxMxL matrix


V1 = permute(ES,[1 3 2]);   % normalize to spherical spreading losses!
V1 = repmat(V1,[1 M 1]);

V2 = permute(B,[3 1 2]);
V2 = repmat(V2,[N 1 1]);

% combine results (by summation in dB)
V = V1 + V2;        % Range-dependence + Azimuth-dependence

% plot slices through volumetric data
if PLOTFLAG
    figure;
    [X,Y,Z] = meshgrid(theta,rho,f*1e-3);
    xslice = [0 30];            % azimuth slice
    yslice = [4.5];              % range slice
    zslice = [20 50 80 100];    % frequency slice
    h = slice(X,Y,Z,V,xslice,yslice,zslice); hold on
    %plot3([theta0 theta0], [rho0 rho0], [0 f(end)]*1e-3,'k','linewidth',2); hold off
    %cRange = get(gca,'clim');
    %set(gca,'clim',[cRange(2)-40 cRange(2)])
    set(gca,'ZDir','reverse')
    set(h,'EdgeColor','none')
    shading interp
    xlabel('azimuth (deg)')
    ylabel('range (m)')
    zlabel('frequency (kHz)')
    title(sprintf('Relative Echo Intensity (dB re %g m)',rhoRef))
    colorbar
end

