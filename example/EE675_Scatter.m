% EE 675 | Dielectric Shell - Scattering/ khaldoon Ishmael 
%

% ------------------------------------------------------------------------
%
% This program is designed to calculate the total Electric Field from 
% the Incident and Scatter Fields.
%
%
% 
% ------------------------------------------------------------------------
clear all %#ok<CLALL>
clc
%
% VARIABLES
%
f = 1000*10^(6); % Frequency
w = 2*pi*f; % Angular Frequency (omega)
epr = 4.0; % Dielectric Constant
ep0 = 8.854*10^(-12); % Dielectric Constant in Air
miu0 = 4*pi*10^(-7); % Permeability
k = w*sqrt(miu0*ep0);
N = 100; % # of sections
lambda = (2*pi)/k; % Wavelength (lambda)
a = 0.25*lambda; % Inner Radius
b = 0.30*lambda; % Outer Radius
r = (a+b)/2; % Distance from Center to Middle of Section
Area = (pi*((b^2)-(a^2)))/N; % Area of Section
area = sqrt(Area/pi); % Radius of Section
phi = (2*pi)/N; % Angle of each section
C = zeros(N); % Create C Matrix to Hold Calculated C-Values
Ei = zeros(1,N); % Create Incident Field Matrix
rho = zeros(N);
%
%
%-------------------------------------------------------------------------
for m = 1:N
    for n = 1:N
        rho(m,n) = 2*(a+(b-a)/2)*sin(phi*(m-n)/2); % Set Rho Values for Different Angles
        if m == n
            C(m,n) = 1+(epr-1)*((1i)/2)*(pi*k*area*besselh(1,2,k*area)-2*1i); % C-Matrix----bessel         
        else
            C(m,n) = (1i*pi*k*area/2)*(epr-1)*besselj(1,k*area)*besselh(0,2,k*abs(rho(m,n))); % C-Matrix
        end
        x(1,m) = r*cos(m*phi); % X-Values
        Ei(1,m) = exp(-1i*k*(x(m))); % Incident Field
    end
end
E = Ei/C; % Total Electric Field

%%%%%%%%%echo width 
% x=r*cos(phi);
% y=r*sin(phi);
% b=exp(-1i*2*pi*x)';
% w(n)= 2*pi*3*abs(exp(1i*2*pi*(x.*cos(phi(n))+y.*sin(phi(n))))*epr-1)*E*area*besselj(1,2*k*area)^2/abs(b(n))^2;



X = linspace(0,180,N/2);
plot(X,abs(E(1:N/2)))
title('Electric Field Distribution'); 
xlabel('Degrees');
ylabel('| E |');
