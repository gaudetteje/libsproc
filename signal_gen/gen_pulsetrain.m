function ts = gen_pulsetrain(fs,tLen,pTime,N0,PT)
% GEN_PULSETRAIN  creates a time series waveform using synthetic pulses
%
% Usage:
%   ts = gen_pulsetrain(fs,tLen,pTime,N0,PT)
% 
% Inputs:
%   fs           - sampling rate of generated waveform [Hz]
%   tLen         - total duration of time series [sec]
%   pTime        - start times for each pulse [sec]
%   N0           - broadband noise level [dBA/sqrt(Hz)]
%   PT           - array of structs containing individual pulse definitions
%     .time      - start and stop times for each component [sec]
%     .freq      - start and stop frequencies for each component [Hz]
%     .modFn     - pulse modulation type ['cw','lfm','sfm','hfm']
%     .winFn     - amplitude shading function [window function string or handle]
%     .gain      - amplitude scaling factor for each component
%     .phase     - initial phase for each component [radians]
%
% Outputs:
%   ts           - struct containing time series data
%     .timestamp - date string when structure was created
%     .fs        - sampling rate
%     .snr       - signal-to-noise ratio(s) of pulses
%     .awgn      - noise level [dB/sqrt(Hz)]
%     .pulse     - structure of pulses as generated
%     .time      - column vector of sample time information
%     .data      - column vector of amplitude information
%     .units     - amplitude units of ts.data
%
% To Do:
% - allow option for custom IF/IA functions
%
% Example:
%   fs = 1e6;
%   tLen = 0.1;
%   pTime = [0.01 0.06];
%   N0 = -40;
%
%   P.time = [0 .025]';
%   P.freq = [100e3 50e3]';
%   P.modFn = 'lfm';
%   P.winFn = 'raisedcos';
%   P.gain = 1;
%   P.phase = 0;
%
%   PT(1) = P;
%   P.freq = [50e3 25e3]';
%   PT(2) = P;
%   X = gen_pulsetrain(fs,tLen,pTime,N0,PT);
%


% initialize time series structure
ts.timestamp = datestr(now());
ts.fs = fs;
ts.awgn = N0;
ts.pulse = PT;
ts.time = (0:1/fs:tLen)';
ts.data = zeros(length(ts.time),1);
ts.units = 'arbitrary';


% iterate over each pulse
for p = 1:length(pTime)

    % check for errors
    assert(length(pTime) == length(PT), 'Length of ''pTime'' must equal length of struct, ''PT''')
    assert(min(min(PT(p).time)) >= 0, 'Cannot specify start times less than 0')
    assert(max(max(PT(p).time)) + pTime(p) <= tLen, 'Cannot specify end times greater than ts.time')


    % generate each multicomponent pulse
    x = gen_ifpulse_multi(fs, PT(p));


    % add current pulse to time series signal
    idx = (1:length(x)) + round(pTime(p)*fs);
    ts.data(idx) = ts.data(idx) + x;
end


%% from Sanderson (2003)
% SNR = 10*log10(2E/N0)
% where, E is echo energy flux:  \int{amp^2} [Pa^2]
% and, N0 is noise power / Hz OR mean(amp^2)/BW [Pa^2]

%% from comp.dsp
% SNR = 10*log10(mean(x.^2) / sigma^2)
% where x is the signal, and sigma is the std. dev. of WGN

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate and add noise (if applicable)
% if SNR ~= Inf
%     
%     % calculate signal energy of initial pulse
%     E0 = mean((x.*pStren(1)).^2);
%     
%     % normalize noise power based on measured E0
%     w0 = sqrt(E0 ./ 10.^(SNR./10)) .* randn(size(ts.data));
%     ts.data = ts.data + w0; % add to time series
% 
%     calculate N0 and ENR for each pulse
%     fprintf('Desired noise level: %f\nActual noise level: %f\nDifference: %f\n\n', N0, mean(w0.^2), N0 - mean(w0.^2))
%     fprintf('ENR = %f\n\n', 10*log10(E0/N0))      % Energy-to-noise ratio (dB)
% end

%% generate by absolute noise level (multiply std of noise)
if N0
    w = 10.^(N0/20) .* randn(size(ts.data));
    ts.data = ts.data + w;
end
