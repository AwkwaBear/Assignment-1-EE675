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
lambda = (2*pi*c)/f;     % wavelength lambda in m

%% Variables %%
lengthAntenna = 0.94;
radiusAntenna = 0.0001;

numberOfSegments = 60;



%% Compute Wavelength %%
lengthLamb = lengthAntenna*lambda;
radiusLamb = radiusAntenna*lambda;



