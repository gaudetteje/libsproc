
% implementation of equations from Urick - Table 3.2 - pg. 43

% piston transducer beam pattern

F = (20e3:10e3:100e3).';
c = 344;
lambda = c./F;
theta = (-90:90);



% piston transducer beam (diameter D in infinite baffle)
D = 0.0047*2;
beta = sin(theta*pi/180);
B1 = (2 * besselj(1,(pi*D./lambda) * beta) ./ (pi*D./lambda * beta)).^2;

figure
plot(theta,db(real(B1),'power'))
axis([-90 90 -40 0])
grid


% 
% % line array of n elements (spacing of d)
% n = 2;
% d = .014;
% B2 = (sin(n*pi*d*sin((theta*pi/180) ./lambda)) ./ ...
%       (n*sin((pi*d./lambda)*sin(theta*pi/180))) ).^2;
% 
% figure
% plot(theta,db(B2))
% grid
