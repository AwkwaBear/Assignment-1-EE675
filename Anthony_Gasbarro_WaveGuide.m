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

analyticalCutoffFreq = (speedOfLight/(2*pi))*(sqrt(((i*pi/(width*h))^2)+((j*pi/(height*h))^2)));

%% Generate A Matrix
Amatrix = -4*eye(nodesTotal,nodesTotal);
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

%% Compute Eigenvalues and Extract
[v, c] = eig(Amatrix);
eigenValue = diag(c);

%% Compute Cutoff Frequency
cutoffFreq = (1/(2*pi))*(speedOfLight/h)*sqrt(-eigenValue);
difference = abs(analyticalCutoffFreq-cutoffFreq);
min = find(abs(difference-min(difference))<0.001);
closeCutoffFreq = cutoffFreq(min(1,:),1);

%% Calculate Wave Velocity
Beta = (2*pi*analyticalCutoffFreq)*(1/speedOfLight)*sqrt(1-(closeCutoffFreq/analyticalCutoffFreq));
waveVelocity = speedOfLight/sqrt(1-(cutoffFreq/analyticalCutoffFreq).^2);

%% Compute Coefficient and input into A-Matrix
coefficient = ((2*pi*closeCutoffFreq)^2)*mu0*E0;
Amatrix = ((coefficient*h^2)-4)*eye(nodesTotal, nodesTotal)+diag1+diag2;


%% Generate output and Graph
angle = 120*pi*sqrt(1-(closeCutoffFreq/analyticalCutoffFreq)^2);
Phi = v(:,min(1,:));
Temp = reshape(Phi, nodesHoriz, nodesVert);
% add two to each dimension for top, bottom, left,and right sides
outputMatrix = zeros(nodesVert+2, nodesHoriz+2);
x = 0:nodesVert-1;
y = 0:nodesHoriz-1;

% place potentials at each border
for i = 1:nodesVert
    outputMatrix(1,1:nodesHoriz+2) = top;
    outputMatrix(i+2,1:nodesHoriz+2) = bottom;
    outputMatrix(i+1, nodesHoriz+2) = right;
    outputMatrix(i+1, 1) = left; 
end
outputMatrix(x+2,y+2) = Temp.';
%fprintf('cutoff frequency=%d\n',cutoffFreq);
fprintf('Ananlytical CutFreq=%d\n',analyticalCutoffFreq);
fprintf('Closest Freq to Ananlytical CutFreq=%d\n',closeCutoffFreq);

imagesc(Temp.');
contourf(Temp.'/angle);

colorbar