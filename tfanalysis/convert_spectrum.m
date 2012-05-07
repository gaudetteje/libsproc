function fd = convert_spectrum(fd,units)
% CONVERT_SPECTRUM  converts the units of a frequency domain structure
%
% FD = convert_spectrum(FD,'units',UNITS) returns a new frequency domain structure
%     with the appropriate units set
% FD = convert_spectrum(FD,'band','full') returns a new frequency domain
%     structure with a (scaled) double sided spectrum
% FD = convert_spectrum(FD,'band','half') returns a new frequency domain
%     structure with a (scaled) single sided spectrum
%
% Available units (case insensitive):
%   (useful for tonal analysis)
%     'Vpk'         - unnormalized peak amplitude [default units for FFT]
%     'Vrms'        - unnormalized rms amplitude
%     'Vpk^2'       - unnormalized peak power
%     'Vrms^2'      - unnormalized rms power
%   (useful for noise floor analysis)
%     'Vpk/rtHz'    - normalized peak amplitude
%     'Vrms/rtHz'   - normalized rms amplitude
%     'Vpk^2/Hz'    - normalized peak power
%     'Vrms^2/Hz'   - normalized rms power
%
% Note:  unnormalized units are not corrected for FFT bin resolution OR equivalent noise power
%        bandwidth of the window
% Note:  dB is calculated separately - units reflect linear magnitude


%%% TBD - convert between full/half spectrum

if ~isfield(fd,'win')
    warning('CONVERT_SPECTRUM:main','No window information present in frequency domain structure - assuming rectangular window')
    fd.win.name = 'rectwin';
end

if ~isfield(fd.win,'NPBW')
    fd.win.NPBW = calc_noisepbw(fd.win);
end

% ensure we have magnitude spectrum to work with
if ~isfield(fd,'mag')
    if ~isfield(fd,'magdb')
        error('No magnitude information present in frequency domain structure')
    else
        warning('CONVERT_SPECTRUM:main','No linear magnitude information found - converting from fd.magdb')
        fd.mag = 10.^(fd.magdb./20);
    end
end

% calculate bin resolution, if unavailable
if ~isfield(fd,'binres')
    fd.binres = fd.freq(2)-fd.freq(1);
end

% normalize units to Vpk (standard FFT output)
switch lower(fd.units)
    case 'vpk'
        mag = fd.mag;
    case 'vrms'
        mag = fd.mag .* sqrt(2);
    case 'vpk^2'
        mag = sqrt(fd.mag);
    case 'vrms^2'
        mag = sqrt(fd.mag .* 2);
    case 'vpk/rthz'
        mag = fd.mag .* sqrt(fd.binres * fd.win.NPBW);
    case 'vrms/rthz'
        mag = fd.mag .* sqrt(fd.binres * fd.win.NPBW * 2);
    case 'vpk^2/hz'
        mag = sqrt(fd.mag .* (fd.binres * fd.win.NPBW));
    case 'vrms^2/hz'
        mag = sqrt(fd.mag .* (fd.binres * fd.win.NPBW * 2));
    otherwise
        error('Unrecognized units in frequency domain structure')
end

% convert from vpk to desired units
switch lower(units)
    case 'vpk'
        fd.mag = mag;
    case 'vrms'
        fd.mag = mag ./ sqrt(2);
    case 'vpk^2'
        fd.mag = mag.^2;
    case 'vrms^2'
        fd.mag = (mag ./ sqrt(2)).^2;
    case 'vpk/rthz'
        fd.mag = mag ./ sqrt(fd.binres * fd.win.NPBW);
    case 'vrms/rthz'
        fd.mag = mag ./ sqrt(fd.binres * fd.win.NPBW * 2);
    case 'vpk^2/hz'
        fd.mag = (mag ./ sqrt(fd.binres * fd.win.NPBW)).^2;
    case 'vrms^2/hz'
        fd.mag = (mag ./ sqrt(fd.binres * fd.win.NPBW * 2)).^2;
    otherwise
        error('Unrecognized units entered')
end

fd.magdb = 20.*log10(fd.mag);
fd.units = units;
