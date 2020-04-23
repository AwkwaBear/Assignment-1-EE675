%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Scattering                                    %%
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
epr = 4.0;
frequency = 1000*10^6;
w = 2*pi*frequency;
k = w*sqrt(mu0*E0);
Lambda = (2*pi)/k;

%% Input Variables %%

numberOfSections = 100;


%% Computed Variables %%

innerRadius = Lambda*0.25;
outerRadius = Lambda*0.3;
r = (innerRadius+outerRadius)/2;
Area = (pi*((outerRadius^2)-(innerRadius^2)))/numberOfSections;
radius = sqrt(Area/pi);
phi = (2*pi)/numberOfSections;

%% Generate Matrices %%

EiMatrix = zeros(1, numberOfSections);
rhoMatrix = zeros(numberOfSections);
CMatrix = zeros(numberOfSections);

for a = 1:numberOfSections
    for b = 1:numberOfSections
        rhoMatrix(a,b) = 2*(innerRadius+(outerRadius-innerRadius)/2)*sin(phi*(a-b)/2); % Set Rho Values for Different Angles
        if a == b
            CMatrix(a,b) = 1+(epr-1)*((1i)/2)*(pi*k*radius*besselh(1,2,k*radius)-2*1i); % C-Matrix----bessel         
        else
            CMatrix(a,b) = (1i*pi*k*radius/2)*(epr-1)*besselj(1,k*radius)*besselh(0,2,k*abs(rhoMatrix(a,b))); % C-Matrix
        end
        x(1,a) = r*cos(a*phi); % X-Values
        EiMatrix(1,a) = exp(-1i*k*(x(a))); % Incident Field
    end
end

%% Compute Total Electric Field %%
totalElectricField = EiMatrix/CMatrix;


%% Graph Output %%
X = linspace(0,180,numberOfSections/2);
plot(X,abs(totalElectricField(1:numberOfSections/2)))
title('Electric Field Distribution'); 
xlabel('Degrees');
ylabel('| E |');
