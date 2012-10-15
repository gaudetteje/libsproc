% Plot sector scan of relative echo intensity, time delay, and/or sensitivity

% Basic assumptions:
%   - Transmission losses are range and frequency dependent and accounted
%     for by only spherical spreading and atmospheric absorption
%   - Transducer is a Piston-type of 4.7 cm radius, as found by Hartley and
%     Suthers (1989)
%   - Receiver(s) are omnidirectional and do not contribute to angular
%     directivity
%   - Echo strength is uniform across all frequencies

%close all
clear
clc

PLOTFLAG = false;

set(0,'DefaultFigureColor','w')
set(0,'DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18)


% region of focus and sensitivity
rho0 = 4.5;               % central region of focus in range [m]
theta0 = 30;             % central region of focus in azimuth [m]
alt_const = -16e-6;     % amplitude latency trading constant [us/dB]
tau = 1e-5;             % delay sensitivity for contour plot

% physical parameters
D = 0.0094;             % transmit aperture
d = 0.014;              % distance between receive elements

% environmental parameters
T = 25;                 % temperature [deg. C]
rh = 50;                % relative humidity [%]
p = 101.325;            % barometric pressure (1 atm)

% plotting parameters
fNum = 401;
rNum = 401;
tNum = 401;
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
TL = TL_abs;%%%%%% + TL_sph;   % TL [dB] = dB + dB

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
    xslice = theta0;            % azimuth slice
    yslice = rho0;              % range slice
    zslice = [20 50 80 100];    % frequency slice
    h = slice(X,Y,Z,V,xslice,yslice,zslice); hold on
    plot3([theta0 theta0], [rho0 rho0], [0 f(end)]*1e-3,'k','linewidth',2); hold off
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


%% compute the effective transfer function at specified focal point in the range/azimuth plane

% find TF at focus point
rIdx = find(rho0 <= rho,1);
thetaIdx = find(theta0 <= theta,1);

if (isempty(rIdx) || isempty(thetaIdx))
    error('Select a focus point within the plotted region')
end
tf = V(rIdx,thetaIdx,:);

% subtract TF from all other points in range-azimuth plane
W = V - repmat(tf,[N M 1]);


% plot slices through volumetric data
if PLOTFLAG
    figure;
    xslice = theta0;        % azimuth slice
    yslice = rho0;        % range slice
    zslice = [20 50 80 100];      % frequency slice
    h = slice(X,Y,Z,W,xslice,yslice,zslice);
    set(gca,'clim',[-40 40])
    set(gca,'ZDir','reverse')
    set(h,'EdgeColor','none')  %'FaceColor','interp')
    %shading interp
    xlabel('azimuth (deg)')
    ylabel('range (m)')
    zlabel('frequency (kHz)')
    title(sprintf('Echo amplitude difference (focal point @ %g m, %g deg)',rho0,theta0))
    colormap(hotcold)
    colorbar
end


%% assign deviation based on error criterion (linear)
E_dbamp = sum(abs(W),3)./L;        % spectrogram correlation deviation (simple summation across frequency)
E_sc = alt_const * E_dbamp;        % convert from dB amplitude error to delay error

% plot normalized L1 norm error surface [dB]
if 1 %PLOTFLAG
    
    dispmeth = 3;       % linear/logarithmic & amplitude (dB) vs. time (microsec)
    plotmeth = 2;       % plot type
    
    figure;
    
    % set intensity scale and depth to match physical parameters
    switch dispmeth
        case 1
            % linear amplitude
            EE = 10.^(-E_dbamp/10);
            cMap = jet;
            cRange = [0 max(max(EE))];
            aRatio = [1 1 1/5];
            units = 'normalized amplitude';
        case 2
            % logarithmic amplitude [dB deviation]
            EE = E_dbamp;
            cMap = flipud(jet);
            cRange = [0 3];
            aRatio = [1 1 1/5];
            units = 'dB';
        case 3
            % linear time scale [microseconds deviation]
            EE = -E_sc*1e6;
            cMap = flipud(jet);
            cRange = [0 40];
            aRatio = [1 1 50];
            units = '\mus';
    end
    
    switch plotmeth
        case 1
            surfc(theta,rho,flipud(EE));
            shading interp
            colorbar
            view(2)
        case 2
            polar3d(flipud(EE),-pi/2,pi/2,rho(1),rho(end),1,'surf'); hold on;
            polar3d(nan(size(EE)),-pi/2,pi/2,rho(1),rho(end)+.1,1,'meshl'); hold off;
            set(gca,'FontSize',18)
            set(gca,'DataAspectRatio',aRatio)
            set(gca,'ylim',rho(end) * [-1.1 1.1])
            set(gca,'xlim',rho(end) * [-0.1 1.1])
            shading interp
            set(gca,'clim',cRange)
            zRange = get(gca,'zlim');
            set(gca,'zlim',[zRange(1)-10 zRange(2)])
            colormap(cMap)
            colorbar('location','westoutside')
            view(2)
        case 3
            polar3d(EE,-pi/2,pi/2,rho(1),rho(end),1,'contour');
    end
    title(sprintf('Minimum L1 error [%s] (%g m, %g deg)',units,rho0,theta0))
end



%% compute TDOA based on 2 spaced out receive elements
rcv = [-d/2 d/2];   % assign coordinates of receive elements
[y0,x0] = pol2cart(theta0*pi/180,rho0);

% calculate TDOA from focal point [s]
tdoa0 = calcArrivalTime([rcv; 0 0], [x0; y0], c);
deltau0 = diff(tdoa0);

% calculate TDOA from all points [s]
[THETA,RHO] = meshgrid(theta',rho);
[Y0,X0] = pol2cart(THETA*pi/180,RHO);
XX = reshape(X0,1,numel(X0));
YY = reshape(Y0,1,numel(Y0));
TDOA = calcArrivalTime([rcv; 0 0],[XX;YY],c);

% compute difference in TDOA from focal point
delTau = diff(TDOA,1,1);
delTau = reshape(delTau,N,M);
ZZ = (delTau - deltau0)*1e6;       % relative difference in time from expected (micro sec)

% plot sector plot showing time difference (in micro sec)
if 1 %PLOTFLAG
    aRatio = [1 1 50];
    figure
    polar3d(flipud(ZZ),-pi/2,pi/2,rho(1),rho(end),1,'surf'); hold on;
    polar3d(nan(size(ZZ)),-pi/2,pi/2,rho(1),rho(end)+.1,1,'meshl'); hold off;
    set(gca,'FontSize',18)
    set(gca,'DataAspectRatio',aRatio)
    set(gca,'ylim',rho(end) * [-1.1 1.1])
    set(gca,'xlim',rho(end) * [-0.1 1.1])
    shading interp
    set(gca,'clim',[-40 40])
    zRange = get(gca,'zlim');
    set(gca,'zlim',[zRange(1)-10 zRange(2)])
    colormap(hotcold)
    colorbar('location','westoutside')
    view(2)
    
    title(sprintf('Time difference of arrival [%s] (%g m, %g deg)','\mus',rho0,theta0))
end


%% combine SC deviation and TDOA deviation

QQ = abs(ZZ) + EE;

if 1 %PLOTFLAG
    aRatio = [1 1 50];
    figure
    polar3d(flipud(QQ),-pi/2,pi/2,rho(1),rho(end),1,'surf'); hold on;
    polar3d(nan(size(QQ)),-pi/2,pi/2,rho(1),rho(end)+.1,1,'meshl'); hold off;
    set(gca,'FontSize',18)
    set(gca,'DataAspectRatio',aRatio)
    set(gca,'ylim',rho(end) * [-1.1 1.1])
    set(gca,'xlim',rho(end) * [-0.1 1.1])
    shading interp
    set(gca,'clim',[0 30])
    zRange = get(gca,'zlim');
    set(gca,'zlim',[zRange(1)-10 zRange(2)])
    colormap(flipud(jet))
    colorbar('location','westoutside')
    view(2)
    
    %title(sprintf('Time difference of arrival [%s] (%g m, %g deg)','\mus',rho0,theta0))
end