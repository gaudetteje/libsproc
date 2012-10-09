close all
clear
clc

% plot sector scan of relative echo intensity, time delay, and/or sensitivity

% user parameters
f = 10.^(4:.1:5);       % frequency points - log scale [Hz]
r = linspace(0.01,10,101);         % range points - linear scale [m]
theta = -90:1:90;       % azimuth points - linear scale [degrees]

% select region of focus
r0 = 1;               % central region of focus in range [m]
theta0 = 0;             % central region of focus in azimuth [m]

% define region of delay sensitivity
tau = 1e-5;             % delay sensitivity for contour plot
alt_const = -16e-6;     % amplitude latency trading constant [us/dB]

%% compute range dependent losses

% calculate frequency dependent absorption and spreading losses
alpha = calcAbsorptionCoef(f,50,25,101.325);
TL_abs = alpha * r;     % TL [dB] => [dB/m * m]

% calculate range dependent spherical spreading losses
TL_sph = ones(numel(f),1) * db(r,'power');



% find total combined 2-way transmission loss vs. range and frequency
TL = TL_abs + TL_sph;   % TL [dB] = dB + dB
E_r = 2 * TL;    % echo strength relative to source level [dB SPL re 1m]

figure(1);
semilogx(r,TL_abs,r,TL_sph)

figure(2);
plot(r,TL)


%% compute beam pattern directivity

