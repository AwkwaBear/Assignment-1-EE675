clearvars;
close all;
clc;

%Computation for Transmission Line Geometry #3
%By: Anthony Gasbarro

%Dimensions of TL
    height = 2; 
    width= 7;

%size of mesh
    h = 0.1;

%Conductor Parameters
phiConductor = 10.00; %potential of the conductor
widthConductor = 1.00;
thicknessConductor = 0; %Assume thickness to be zero

heightDielectric = 1.00;
er = 9.6;

%Contour Distance from Center Conductor
contourdh = 5;
contourdv = 7;
rmax = (height/h)-1;
cmax = (width/h)-1;
nodesTotal = rmax*cmax;


%generate A Matrix with diagonal 1's and -4's
Amatrix = zeros(nodesTotal, nodesTotal); 

for x = 1:rmax
    
end