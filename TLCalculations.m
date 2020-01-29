%Computation for Transmission Line Geometry #3
%By: Anthony Gasbarro

%Mesh Dimensions
    heightMesh = 2; 
    widthMesh = 7;

%size of mesh
    size = .1;

%Conductor Parameters
voltageConductor = 10.00;
widthConductor = 1.00;
thicknessConductor = 0; %Assume thickness to be zero

heightDielectric = 1.00;
er = 9.6;

%Contour Distance from Center Conductor
contourdh = 5;
contourdv = 7;
rmax = heightMesh/(size-1);
cmax = widthMesh/(size -1);
nodesTotal = rmax*cmax;