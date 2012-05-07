function rmsnoise = total_spectral_noise(fd,varargin)
% TOTAL_SPECTRAL_NOISE  calculates the total RMS noise from a frequency
% domain structure
%
% RMSNOISE = calc_spectral_noise(FD) returns a scalar value of the noise
%     contained within the entire spectrum
% RMSNOISE = calc_spectral_noise(FD,FRANGE) returns a scalar value of the
%     noise contained within a band limited spectrum - FRANGE = [fmin fmax]
%

% convert to power spectrum and use normalized bin widths
fd = convert_spectrum(fd,'Vrms^2/Hz');

% select frequency range (if necessary)
if nargin > 1
    mag = [];
    for n = 1:nargin-1
        frange = varargin{n};
        mag = [mag fd.mag(fd.freq >= frange(1) & fd.freq <= frange(2))];
    end
%%%%%%%% TBD need to warn if points are missing!
else
    mag = fd.mag;
end

% integrate the points in each 1Hz bin
rmsnoise = sqrt(sum(mag).*(2*fd.binres));
