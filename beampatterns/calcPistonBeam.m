function varargout = calcPistonBeam(D,f,varargin)
% calcPistonBeam  Computes the beam pattern response for a standard piston
% model of fixed aperture in an infinite baffle
%
% B = calcPistonBeam(D,f) returns the beam pattern vs. 
%
% Default parameters:
%   c = 344             speed of sound in air [m/s]
%   theta = (-90:2:90)  azimuth [degrees]
%
% Example 1 - Calculate the 1kHz beam pattern for an aperture of a 2 cm diameter
%   >> D = 0.02;
%   >> f = 1e3;
%   >> B = calcPistonBeam(D,f);
%
% Example 2 - Calculate the beam pattern across two decades of frequency.
%   Use specific angles and sound speed for 15 degrees Celcius. 
%   >> D = 0.02;
%   >> f = 10.^(2:0.1:4)
%   >> theta = (-45:5:45);
%   >> c = calcSoundSpeed(15);
%   >> B = calcPistonBeam(D,f,theta,c);
%
% Note:  The beam pattern is symmetrical about 0 degrees in both azimuth
% and elevation.
%
% See following references for more information:
%   Urick - Table 3.2 - pg. 43

% assign default values
theta = (-90:1:90).';     % angle from endfire [degrees]
c = 344;                % speed of sound [m/s]

switch nargin
    case 4
        theta = varargin{1};
        c = varargin{2};
    case 3
        theta = varargin{1};
    case 2
    otherwise
        error('Incorrect number of parameters entered')
end

% force angles into column vector
theta = theta(:);

% convert frequency to wavelength
lambda = c./f;

% compute theoretical piston beam pattern
beta = sin(theta*pi/180);
B = (2 * besselj(1,beta * (pi*D./lambda)) ./ (beta * (pi*D./lambda))).^2;

% generate plot if no output arguments present
switch (nargout)
    case 0
        figure
        plot(theta,db(real(B),'power'))
        axis([theta(1) theta(end) db(min(min(B)),'power') db(max(max(B)),'power')])
        grid on
    case 1
        varargout{1} = B;
    case 2
        varargout{1} = B;
        varargout{2} = beta;
    otherwise
        error('Incorrect number of output parameters')
end
