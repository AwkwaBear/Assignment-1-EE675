%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% TM Waveguide                                  %%
%% Anthony Gasbarro                              %%
%% EE 675                                        %%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clearvars;
close all;
clc;


%% Constants %%

speedOfLight = 2.997925 *10^8; %Velocity of Light through air
E0 = 8.854*10^(-12); %permittivity of free space (F/m)
mu0 = 4*pi*10^(-7);     %permeability of free space in N/(A^2)

%% Adjustable parameters %%

h = 0.5; %resolution
width = 10/h; %Width of area
height = 20/h; %Height of area
%set potentials for 
top = 0;
bottom = 0;
right = 0;
left = 0;

%% Compute Grid Dimensions %%
nodesVert = height-1;
nodesHoriz = width-1;
nodesTotal = nodesVert * nodesHoriz;
%set node iterators to node (1,1)
i = 1; 
j = 1;

