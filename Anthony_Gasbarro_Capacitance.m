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
diag = sqrt(DA)/(2*E0*sqrt(pi));
hDist = Distance/h;


%% Matrix values %%
count = 1;
z = 1;
stop = N;
rBottom = zeros(n);
rTop = zeros(n);
diagonal = sqrt(DA)/(2*E0*sqrt(pi));
seta = zeros(2*N, 2*N);
q = zeros(1, hDist);
distanceVector = zeros(2*N);

%% Create Distance Matrix %% 

for k = 1:hDist
    for x = 1:n
        for y = 1:n
            for i = 1:n
                for j = 1:n
                    rTop(i,j) = sqrt(((abs(y-i))*Length)^2 + ((abs(x-j))*Length)^2);
                    rBottom(i,j) = sqrt((z*h)^2 + (rTop(i,j))^2);
                end
            end
            RTOP = rTop(:)';
            RBOTTOM = rBottom(:)';
            row1 = horzcat(RTOP, RBOTTOM);
            row2 = horzcat(RBOTTOM, RTOP);

            if count <= N
                distanceVector(count, :) = row1;
                distanceVector(N+count, :) = row2;
                count = count+1;
            end
        end
    end
    if z < hDist
        count = 1;
    end
    Vector = distanceVector + eye(2*N);
    B = 1./Vector;
    A = (ones(2*N)*coefficient) - eye(2*N)*coefficient; 
    A1 = eye(2*N)*diagonal; 
    A2 = A + A1; 
    Amatrix = A2.*B;
    pos = ones(N, 1);
    neg = -ones(N, 1);
    Bmatrix = vertcat(pos, neg);
    rho = Amatrix\Bmatrix;
    halfpos = rho(1:N);
    halfneg = rho(N+1:2*N); % Calculating total rho Values from -
    R1 = reshape(halfpos,n,n);
    R0 = -R1;
    q(1,z) = sum(halfpos*SubArea);
    z = z+1;
end

%% Compute Capacitance Ratio %%
dist = h:h:Distance;
dwratio = dist./Width;
prod = 1./dwratio.*Width;
C0 = prod.*Width.*E0; % analytical capacitance
cap = q./(Vdiff-(-Vdiff));

capratio = cap./C0;


%% Graph Final Outputs %%

%display 2D charge density plot where left side is top plate, right side
    %is bottom plate, they are folded open from each other over center line
    %left to right is length, top to bottom is width

figure(1)
subplot(2,1,1)
imagesc(R1)
colorbar
hold on
subplot(2,1,2)
imagesc(R0)
colorbar
hold off

%display 3D charge density plot
figure(2)
surf(R1)
hold on
surf(R0)
hold off
colorbar

%plot normalized capacitance C/C0 for d/W
figure(3)
plot(dwratio,capratio,'-o')
 title('Normalized Capacitance vs d/W for Delta Function Basis')
 xlabel('d/W')
 ylabel('Normalized Capacitance C/C0')
    
% title('Parallel Plate Capacitance with The Increase in Distance')
% xlabel('d/W')
% ylabel('C/C0')
