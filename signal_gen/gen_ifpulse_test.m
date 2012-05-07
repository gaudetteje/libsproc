% This script demonstrates usage examples for gen_ifpulse()

clc
%clear
%close all


% user defined parameters
fs = 500e3;%            % sampling rate of generated waveform [Hz]
T = 0.005;              % pulse length [sec]
f0 = 110e3;             % initial frequency [Hz]
f1 = 10e3;              % final frequency [Hz]
phi0 = -pi/2;           % inital phase [radians]
modfn = 'lfm';          % predetermined frequency modulation type
ampfn = 'rcos';         % windowed amplitude function


% generate time series vector
t = (0:1/fs:T)';
B = abs(f1-f0);


% define instantaneous amplitude vector
A = ones(length(t),1);
switch ampfn
    case 'rcos'
        idx = round(0.0625*length(A));
        A(1:idx+1) = A(1:idx+1) .* cos(pi/2*(0:1/idx:1) - pi/2)';
        A(end-idx:end) = A(end-idx:end) .* cos(pi/2*(0:1/idx:1))';
    case 'none'
end


% define instantaneous frequency, IF(t), and its integral, phi_ref(t)
switch modfn
    case 'cw'
        %%% CW
        IF = f0*ones(size(t));
        phiref = f0.*t;
    
    case 'lfm'
        %%% LFM
        IF = f0*ones(size(t)) + (f1-f0)/(T).*t;
        phiref = f0.*t + 0.5*(f1-f0)/(T).*t.^2;
    
    case 'sfm'
        %%% SFM
        fm = 5000;
        IF = B./2*sin(2*pi*fm*t) + f0+B/2;
        phiref = -B/(4*pi*fm).*cos(2*pi*fm*t) + f0+B/2*t;  % wrong phase
    
    case 'hfm'      % only valid for downward sweep
        %%% HFM
        a = T*(f0*f1)/B;
        b = T*f1/B;
        IF = a./(t+b);
        phiref = a*log(t+b); % + 0.775/2;    % for some reason, phase is off from 0 radians - the amount depends on T and fs
end


% generate phase modulated pulse
[x,phi] = gen_ifpulse(fs,IF,phi0,A);

% generate reference waveform
xref = A .* exp(1i*2*pi*phiref + 1i*phi0);


% plot time series
if 1
    figure;
    subplot(2,1,1)
    plot(t,real(x),'b',t,real(xref),'r');
    grid on;
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    title('x(t)')
    legend('generated x(t)','reference x(t)')
    axis([t(1) t(end) -1.1*max(A) 1.1*max(A)])
    
    subplot(2,1,2)
    plot(t,real(xref)-real(x))
    grid on;
    xlabel('Time (seconds)')
    ylabel('Amplitude')
    title('Integration Error:  x_{ref}(t) - x(t)')
end


% check that the IF was defined properly
if 1
    figure;
%    subplot(2,1,1)
    plot(t,IF,'g',t(1:end-1),fs*diff(phi),'--b')%,t(1:end-1),fs*diff(phiref),'--r')
    grid on;
    xlabel('Time (seconds)')
    ylabel('Instantaneous Frequency (Hz)')
    title('d\phi / dt')
    legend('Desired d\phi / dt','Approximated d\phi / dt', 'Reference d\phi_{ref} / dt')
    axis([0 T 0 fs/2])
    
%     subplot(2,1,2)
%     plot(t(2:end),IF(2:end)-fs*diff(phi),'r')
%     title('Integration Error:  \phi_{ref} - \phi')
%     grid on;
%     xlabel('Time (seconds)')
%     ylabel('Frequency (Hz)')
end


% perform TF analyses of waveform
if 1
    % spectral analysis parameters
    nfft = 112;
    overlap = 100;
    wind = window(@blackman7,nfft);
    
    % plot the periodogram
    figure;
    pwelch(x,wind,overlap,nfft,fs);
    
    % plot the spectrogram
    figure;
    spectrogram(real(x),wind,overlap,nfft,fs,'yaxis')
    
    fprintf('Calculated power and energy in pulse:  %f, %f\n', mean(x.^2), sum(x.^2));

end

