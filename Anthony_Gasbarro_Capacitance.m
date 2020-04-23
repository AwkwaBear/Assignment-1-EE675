%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Capacitance                                   %%
%% Anthony Gasbarro                              %%
%% EE 675                                        %%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clearvars;
close all;
clc;


%% Constants %%

speedOfLight = 2.997925*10^8; %Velocity of Light through air
E0 = 8.854*10^(-12); %permittivity of free space (F/m)
mu0 = 4*pi*10^(-7);     %permeability of free space in N/(A^2)
diag = sqrt(DA)/(2*ep0*sqrt(pi));

%% Variables %%
Vdiff = 1;
Width = 1;
Distance = 1;
n = 15;
h = 0.1;

%% Computed Variables %%
Area = Width^2;
N = n^2;
Length = Width/n;
SubArea = Length^2;
coefficient = SubArea/(4*pi*E0);
DA = pi*(Length/2)^2;
hDist = Distance/h;


%% Matrix values %%
count = 1;
stop = N;
rBottom = zeros(n);
rTop = zeros(n);
diagonal = sqrt(DA)/(2*E0*sqrt(pi));
seta = zeros(2*N, 2*N);
q = zeros(1, hDist);
distanceVector = zeros(2*N);




