function Y = generallinearfilt(data,hmatch,alpha)
% GENERALLINEARFILT  calculates and applies the general linear filter with
% parameter alpha
%
% Y = GENERALLINEARFILT(DATA,REP,ALPHA) transforms the signal, S(f), into
% the generalized linear filter, H(f), with parameter alpha.  Returns
% filtered output from DATA.
%
% H(f) = |S(f)|^(2*alpha-1) exp(-j*angle(S(f)))
%
% When alpha=1, H(f) is the matched filter.  When alpha=0, H(f) is the
% inverse filter.  The GLF is defined for the alpha range 0 to 1.  Stable
% filters will put a minimum bound of ~0.3 on alpha.

% ensure good values for alpha
assert(alpha >= 0, 'alpha parameter is out of bounds:  0 <= alpha <= 1')
assert(alpha <= 1, 'alpha parameter is out of bounds:  0 <= alpha <= 1')
if alpha < 0.3
    warning('GLF:inputerr','alpha should be set between 0.3 and 1.0 for best results')
end

% make the replica a column vector
hmatch = hmatch(:);

% ensure the replica is unity norm
hmatch = hmatch./norm(hmatch);

% calculate GLF impulse response
nfft = length(hmatch);          % set number of points equal to hmatch
S = fft(hmatch,nfft);           % calculate frequency spectrum
H = abs(S) .^ (2*alpha-1) .* exp(-1i*angle(S));
hd = real(ifft(H));             % convert spectral level to impluse response

% estimate TF coeffs
[b,a] = prony(hd,nfft,nfft);

% apply filter to data
Y = filter(b,a,data);
