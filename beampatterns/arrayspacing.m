close all

opts = {'b','r','g','m','k','c'};

%%%%%%%%%
alpha = (1:90);         % maximum incident angle [deg]
D = 0.0254 * [1 0.75 0.5 0.25];       % maximum array element spacing [m]
F = [10e3 20e3 50e3 600e3];          % maximum frequency of interest [Hz]
c = 1500;                        % speed of sound [m/s]

%%%%%%%%
beta = (90-alpha)*pi/180;       % find opposite angle and convert to radians
DF = c/2 ./ cos(beta);                 % calculate upper bound on min element spacing x max frequency



%%%%%%%%
% show D given F
figure;
hold on
for n=1:length(F)
    plot(alpha,1e2.*DF ./ F(n),opts{n})
end
%axis([0 90 0 5])
grid on
legend('10kHz','20kHz','50kHz','100kHz','location','northeast')
title('Maximum Element Spacing, D, given F')
xlabel('Maximum Incident Angle (deg)')
ylabel('Maximum Element Spacing (cm)')

% show F given D
figure;
hold on
for n=1:length(D)
    plot(alpha,1e-3.*DF ./ D(n),opts{n})
end
%axis([0 90 0 120])
grid on
legend('1"','0.75"','0.5"','0.25"','location','northeast')
title('Maximum Frequency, F, given D')
xlabel('Maximum Incident Angle (deg)')
ylabel('Maximum Frequency (kHz)')

