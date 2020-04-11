%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% TM Waveguide                                  %%
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

%% Adjustable parameters %%

h = 0.5; %resolution
width = 10/h; %Width of area
height = 20/h; %Height of area

% set potentials for sides
top = 0;
bottom = 0;
right = 0;
left = 0;

%% Compute Grid Dimensions %%
nodesVert = height-1;
nodesHoriz = width-1;
nodesTotal = nodesVert * nodesHoriz;

%set node iterators to node (1,1)
i = 1; %horizontal direction iterator
j = 1; %vertical direction iterator

anafreq = (speedOfLight/(2*pi))*(sqrt(((i*pi/(width*h))^2)+((j*pi/(height*h))^2)));

%% Generate A Matrix
Amatrix = -4*eye(nodes,nodes);
% Create Diagnol rows of 1's in A Matrix
for i = 1:nodesHoriz-1
    for j = 1:nodesVert
        diag1(i+(j-1)*nodesHoriz, i+(j-1)* nodesHoriz+1) = 1;
        diag1(i+(j-1)*nodesHoriz+1, i+(j-1)*nodesHoriz) = 1;        
    end
end

for i = 1:nodesHoriz
    for j = 1:nodesVert-1
        diag2(i+(j-1)*nodesHoriz, i+j*nodesHoriz) = 1;
        diag2(i+j*nodesHoriz, i+(j-1)*nodesHoriz) = 1;
    end
end

Amatrix = Amatrix + diag1 + diag2;

[v, c] = eig(Amatrix);
eigenValue = diag(c);

