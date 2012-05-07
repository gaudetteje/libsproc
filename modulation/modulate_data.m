function varargout = modulate_data(X,Fs,Fc,varargin)
% MODULATE_DATA  modulates real or complex basebanded time series
%
% Y = modulate_data(X,Fs,Fc) takes a baseband signal, X, sampled at rate,
%     Fs, and modulates the signal to Fc.  (Fc < Fs/2)
% Y = modulate_data(X,Fs,Fc,I) interpolates X to sampling rate Fs*I before
%     modulating to Fc.  (Fc > Fs/2)
% [Y,I] = modulate_data(...) also returns the interpolation factor, I
% modulate_data(...) will plot the signal before and after modulation
%
% Notes:
%
% X must be a LxN matrix for L samples and N channels of equal length

I = 1;  % default to no interpolation

switch nargin
    case 4
        I = varargin{1};
        if mod(I,1)
            I = round(I);
            warning('Interpolation factor, I, must be a positive integer; Rounding to nearest whole number')
        end
end

L = size(X,1);          % signal length [samples]
C = size(X,2);          % number of channels
if C > L, warning('X has more rows (channels) than (columns) samples.  Proceeding anyway...'), end
if (Fc > I*Fs/2), warning('Fc must be less than the Nyquist rate.  Aliasing is likely to occur.'), end

% Interpolate (zero pad & LP filter)
if I > 1
    X_hat = zeros(L*I,C);
    for ch = 1:C
        X_hat(:,ch) = interp(X(:,ch),I);
    end
else
    X_hat = X;
end

% Create frequency shift vector
t = (0:(I*L)-1).' ./ (I*Fs);
t = repmat(t,1,C);
shifter = exp(1i*2*pi*Fc*t);

% Modulate to analytic signal
Xa = shifter .* X_hat;

% Remove imaginary components to yield 1/2 amplitude real signal
Y = real(Xa);



%
% assign output variables or plot results
%
switch nargout
    case 1
        varargout{1} = Y;
    case 2
        varargout{1} = Y;
        varargout{2} = I;
    case 0      % plots for sanity check and debugging
        wind = @rectwin;
        %wopt = 40;
        
        winX = window(wind,L); %,wopt);
        nfft = 1024;

        Hx = fftshift(fft(winX.*X,nfft))./nfft;

        Fx = Fs.*(0:nfft-1).'./nfft;
        Fx = [-flipud(Fx(1:end/2+1)) ; Fx(2:end/2)];
        Fx = Fx + Fc;

        % Basebanded data spectrum
        figure
        plot(Fx, db(abs(Hx)))
        grid on;
%        set(gca,'XLim',[Fx(1) Fx(end)]);
        title(sprintf('Original (basebanded) signal'))
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
        
        
        
        winY = window(wind,L*I); %,wopt);
        nfft = nfft * I;
        
        Hy = fftshift(fft(winY.*Y,nfft))./nfft;
        
        Fy = Fs.*(0:nfft-1).'./nfft;
        Fy = [-flipud(Fy(1:end/2+1)) ; Fy(2:end/2)];
        Fy = Fy .* I;
        
        % Modulated data spectrum
        figure
        plot(Fy, db(abs(Hy)))
        grid on;
%        set(gca,'XLim',[Fy(1) Fy(end)]);
        title(sprintf('Modulated (carrier frequency) signal'))
        xlabel('Frequency (Hz)')
        ylabel('Magnitude (dB)')
end
