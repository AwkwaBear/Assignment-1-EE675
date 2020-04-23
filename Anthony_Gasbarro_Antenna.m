%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Antenna                                       %%
%% Anthony Gasbarro                              %%
%% EE 675                                        %%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clearvars;
close all;
clc;


%% Constants %%

c = 2.997925*10^8; %Velocity of Light through air
E0 = 8.854*10^(-12); %permittivity of free space (F/m)
mu0 = 4*pi*10^(-7);     %permeability of free space in N/(A^2)
eta0 = 120*pi;           %impedance in free space eta0=377
eps0 = 8.854*10^(-12);   %permittivity of free space in F/m
lambda = (2*pi*c)/(1*10^9);     % wavelength lambda in m
B0 = (2*pi)/lambda;   % phase constant beta0
Vg = 1;                  %voltage difference between gap

%% Variables %%
lengthAntenna = 0.94;
radiusAntenna = 0.0001;

numberOfSegments = 60;

%% Compute Variable Based Values %%
lengthLamb = lengthAntenna*lambda;
radiusLamb = radiusAntenna*lambda;
lengthOfSegment = lengthLamb/numberOfSegments;

segmentCenter = (-lengthLamb/2 + lengthOfSegment/2):lengthOfSegment:(lengthLamb/2 - lengthOfSegment/2);
segmentCenter(numberOfSegments + 1) = lengthLamb/2 - lengthOfSegment;

%% Generate A and B Matrix %%

Amatrix = zeros(numberOfSegments+1, numberOfSegments+1);
Bmatrix = (((Vg/(2*eta0))*1j)*sin(B0*abs(segmentCenter)))';

     %Assign A matrix values
     for i = 1:numberOfSegments+1

         R = sqrt(radiusLamb^2 + (segmentCenter-segmentCenter(i)).^2);
         Amatrix(i, :) = (exp(-1j*B0*R)./R)*lengthOfSegment;
     end
    %Last A matrix Column 
        Amatrix(:,numberOfSegments+1) = cos(B0*segmentCenter)';
        
%% Compute Current %%

I = (Amatrix)\Bmatrix*4*pi;
I = I(1:numberOfSegments);

AntennaLength = segmentCenter(1:numberOfSegments);

Icenter = I(abs(numberOfSegments/2));
Z = Vg/(Icenter);

%% Print Computation Results %%
fprintf('Completed Computation \n')
fprintf('Length: %.5f m \nRadius: %.5f m \nSections: %.8f m\n', lengthLamb, radiusLamb, lengthOfSegment)
fprintf('Gap Potential: %.2f V \nWavelength: %.2f m\n\n', Vg, lambda)

%% Plot final Results %%
plot(abs(I),AntennaLength/lambda)
title('Antenna Current Density Distribution')
ylabel('Section of Antenna-Both Rods(m)') 
xlabel('Current Density (A/(m^2))') 
