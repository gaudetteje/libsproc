function calcHornBeam(D1,D2,L,f,varargin)
% calcHornBeam  Computes the beam pattern response for an acoustic pressure
%    horn given the dimensions (in meters) and frequencies of interest
% 
% B = calcHornBeam(D1,D2,L,F) returns the beam pattern at frequencies, F,
%     for a horn with throat diameter D1, mouth diameter D2, and length, L
% B = calcHornBeam(D1,D2,L,F,theta,beta,c) overrides the default parameters below
%
% Note:  Unless specified otherwise, the horn is conical in shape, normally
%   truncated, and has the following default values.
%
%  theta = -90:1:90;        % angular coverage range [deg]
%  beta = 90;               % angle of truncation [deg]
%  c = 344;                 % speed of sound in the medium [m/s]
%
% See the following references for more information:
%   N.H. Fletcher and S. Thwaites, (1988) "Obliquely truncated simple
%   horns: Idealized models for vertebrate pinnae." Acustica, Vol. 65,
%   pp. 194-204

%TBD

% assign default values
theta = (-90:1:90).';       % angle from endfire [degrees]
c = 344;                    % speed of sound [m/s]
mode = 2;                   % horn type (i.e. parabolic, conical, exponential)

switch nargin
    case 7
        a = varargin{1};
        theta = varargin{2};
        beta = varargin{3};
        c = varargin{4};
    case 6
        a = varargin{1};
        theta = varargin{2};
        beta = varargin{3};
    case 5
        a = varargin{1};
        theta = varargin{2};
    case 4
        a = varargin{1};
    case 3
    otherwise
        error('Incorrect number of parameters entered')
end

assert(0 <= beta && beta <= 90, 'angle of truncation must be between 0 and 90 [deg]!')

% force angles into column vector
theta = theta(:);

% compute wave number, k, and mouth radius, a
k = 2*pi*f/c;
a = D2/2;
alpha = L;      % assumes conical horn (for now)

%% compute the horn impedance parameters
switch (mode)
    % parabolic horn shape
    case 1
        
    % conical horn shape
    case 2
    
    % exponential horn shape
    case 3
        
    otherwise
        error('Unknown mode for horn shape specified')
end

Z12 = ;
Z22 = ;
ZR = ;

%% compute the spherical vs. planar wave correction factor
Falph = sin(0.5*k*a*tan(alpha/2)) ./ (0.5*k*a*tan(alpha/2));

%% compute the magnitude response [dB] along the geometrical axis of the horn
G = 20*log10((Z12 / (Z22 + ZR)) * (1 + tanh(k*a)) * Falph );

%% attenuate response by the corresponding angle (obliquity function)
beta = k*a*sin(theta);
if beta = 90;
    Hbeta = 2 * real(besselj(1,beta)) ./ beta;     % this form reduces to a piston transducer in infinite baffle!
else
    Hbeta = ;
end
Hbeta(isnan(Hbeta)) = 1;    % fix issue with 0/0

%% combine effect of both gain function, G, and angular function, H
B = G + 20*log10(Hbeta);      % summation for case of G and Hbeta specified in dB
