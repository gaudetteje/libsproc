function varargout = baseband_data(X,Fs,Fc,BW,varargin)
% BASEBAND_DATA  applies basebanding to a real or complex time series signal
%
% Y = baseband_data(X,Fs,Fc,BW) takes a modulated signal, X, sampled at rate,
%     Fs, and basebands the signal centered at Fc with bandwidth, BW.
% Y = baseband_data(X,Fs,Fc,BW,D) applies the decimation factor, D.
% Y = baseband_data(X,Fs,Fc,BW,D,N) uses N points for bandpass filter
% [Y,D] = baseband_data(...) also returns the decimation factor, D
%
%
% Notes:
%
% X must be a LxN matrix for L samples and N channels of equal length
% D defaults to an integer yielding at least fs = 2.5 x fmax
% N defaults to 512


% Default parameters
N = 512;    % number of filter coefficients

switch nargin
    case 4
        D = ceil(Fs/(2.5*BW/2));
        fprintf('Decimating basebanded signal by %d\n',D)
    case 5
        D = varargin{1};
    case 6
        D = varargin{1};
        N = varargin{2};
end

L = size(X,1);          % signal length [samples]
C = size(X,2);          % number of channels
M = round(N/2);         % filter delay [samples]
if C > L, warning('X has more rows (channels) than (columns) samples.  Proceeding anyway...'), end



%
% Compute the scaled band limits for the fir1 function
%
ws = [Fc-BW/2 Fc+BW/2];
ws = ws/(Fs/2); % scale the band limits relative to the Nyquist rate


%
% Compute the FIR filter coefficients
%
B = fir1(N,ws);


%
% Apply the bandpass filters to the replica data
%
X_hat = [X ; zeros(M,C)];                   % pad end of signal with zeros
X_hat = filter(B,1,X_hat);                  % filter band of interest
X_hat = X_hat(M+1:end,:);                   % remove the bandpass filter lag
X_ana = hilbert(X_hat);                     % construct analytic signal


%
% Shift spectral content and decimate
%
t = (0:L-1).'/Fs;
t = repmat(t,1,C);
shifter = exp(-1i*2*pi*Fc*t);               % spectral vector to shift data to DC

Y = X_ana .* shifter;                        % shift data to DC
Y = Y(1:D:end,:);                           % decimate by factor D

%
% assign output variables or plot results
%
switch nargout
    case 1
        varargout{1} = Y;
    case 2
        varargout{1} = Y;
        varargout{2} = D;
    case 0      % plots for sanity check and debugging
        wind = @rectwin;
        wopt = 40;
        winX = window(wind,L); %,wopt);
        winY = window(wind,ceil(L/D)); %,wopt);
        nfft = 1024;

        Hx = fftshift(fft(winX.*X,nfft))./nfft;
        Hy = fftshift(fft(winY.*Y,nfft))./nfft;

        Fx = Fs.*(0:nfft-1).'./nfft;
        Fx = [-flipud(Fx(1:end/2+1)) ; Fx(2:end/2)];
        Fy = Fs.*(0:nfft-1).'./nfft;
        Fy = [-flipud(Fy(1:end/2+1)) ; Fy(2:end/2)];
        Fy = Fc + Fy./D;

        % Modulated data spectrum
        figure
        plot(Fx, db(abs(Hx)))
        grid on;
        axis([Fx(1) Fx(end) -100 0])
        title(sprintf('Original (modulated) signal'))
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')

        % Bandpass filter response
        figure
        freqz(B,1,nfft,Fs);
        title(sprintf('Bandpass Filter Frequency Response'))
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')

        % Basebanded data spectrum
        figure
        plot(Fy, db(abs(Hy)))
        grid on;
        axis([Fy(1) Fy(end) -100 0])
        title(sprintf('Baseband (demodulated) signal'))
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
end
