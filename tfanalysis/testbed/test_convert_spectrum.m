clc
clear
close all


% define test signal
ts.fs = 1e6;

A1 = 2*sqrt(2);
A2 = 0.5;
S = 1.66e-15;                   % define noise power spectral density

x1 = A1 * sin(2*pi*2048*(1:ts.fs)*2^-15).';
x2 = A2 * cos(2*pi*4096*(1:ts.fs)*2^-15).';
w = sqrt(S * ts.fs)/1.253 .* randn(ts.fs,1);       % why do we need this fudge factor here? rms(ts.data) shows correct noise level...
%assert(abs(Arms - rms(ts.data)) < 0.01 * Arms, 'Time series is not within tolerance of specified RMS level')   % ensure we have the correct noise power density

% combine tones and noise
ts.data = x1 + x2 + w;

% calculate initial spectrum
fd = calc_spectrum(ts,32768,'rectwin',32);

% loop through all available units
units = {'Vpk','Vrms','Vpk^2','Vrms^2','Vpk/rtHz','Vrms/rtHz','Vpk^2/Hz','Vrms^2/Hz'};

for n=1:length(units)
    % convert to desired units
    fd = convert_spectrum(fd,units{n});

    % check tone level
    fprintf('\nTone 1 is %g %s (%g %s)\n',fd.mag(2049),fd.units,fd.magdb(2049),['dB' fd.units])
    fprintf('Tone 2 is %g %s (%g %s)\n',fd.mag(4097),fd.units,fd.magdb(4097),['dB' fd.units])
    
    % check noise level
    fprintf('Noise spectral density is %g %s (%g %s)\n',mean(fd.mag),fd.units,mean(fd.magdb),['dB' fd.units])%total_spectral_noise(fd))
    fprintf('Total RMS noise %g Vrms\n',total_spectral_noise(fd));

    % plot magnitude spectrum
    figure
    plot(fd.freq,fd.magdb);
    ylabel(sprintf('Magnitude (%s)',fd.units))
    xlabel('Frequency (Hz)');
    grid on;
    hold on;

end
tilefigs(2,4)

% convert back to initial units to make sure no errors occurred
fd = convert_spectrum(fd,units{1});

% recheck tone level
fprintf('\nTone 1 is %g %s (%g %s)\n',fd.mag(2049),fd.units,fd.magdb(2049),['dB' fd.units])
fprintf('Tone 2 is %g %s (%g %s)\n',fd.mag(4097),fd.units,fd.magdb(4097),['dB' fd.units])

% recheck noise level
fprintf('Noise spectral density is %g %s (%g %s)\n',mean(fd.mag),fd.units,mean(fd.magdb),['dB' fd.units])%total_spectral_noise(fd))
fprintf('Total RMS noise %g Vrms\n',total_spectral_noise(fd));

