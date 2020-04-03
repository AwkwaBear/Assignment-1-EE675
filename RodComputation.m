%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computation for ROD                           %%
%% Anthony Gasbarro                              %%
%% EE 675                                        %%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

clearvars;
close all;
clc;

%% Poll User Input %%

prompt = 'How many subdivisions does the rod have? \n'; %prompt user for input
subdivisions = input(prompt); %gather input to create number of rod subdivisions

%% Adjustable parameters %%

length = 1; %length of rod in meters
radius = 0.001; %radius of rod in meters
potential = 1; %Potential at center of section in Volts

%% Constant Parameters %%
E0 = 8.854*10^(-12); %permittivity of free space (F/m)
pi = 3.14159;            %value of pi
mu0 = 4*pi*10^(-7);     %permeability of free space in N/(A^2)

%% Compute section lengths %%

h = 1/subdivisions;
sections = length/h+1;

%% Generate A Matrix %%

%This section

Amatrix = zeros(sections);
for i = 1:sections
        for j = 1:sections
                Amatrix(i,j) = (2*pi*radius*h)/abs(i-j);
        end
end
for i = 1:sections
    Amatrix(i,i) = (4*pi*radius*log(h/radius));
end

%% Generate B Matrix %%

Bmatrix = ones(sections, 1)*4*pi*E0;

%% Calculate C Matrix %%

Cmatrix = Amatrix\Bmatrix;

%% Calculate Charge %%

charge = Cmatrix*2*pi*radius*length;

%% Determine Axis and Plot %%

y = [charge;charge(numel(charge))];
x = (0:1/(numel(charge)):1);

stairs(x,y);
title('Rod Charges');
ylabel('Length in meters');
xlabel('charges [l/\rho]');

