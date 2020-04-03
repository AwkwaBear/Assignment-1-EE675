% EE 471 | Waveguide
%
% Khaldoon Ismaeel 
%
% ------------------------------------------------------------------------
%
% This program is designed to calculate the cutoff frequency, velocity, and
% show the E-Field.
%
%
% 
% ------------------------------------------------------------------------
clear
clc
% prompt = {'What is the value of h? ', 'Enter Width ', 'Enter Height '}; % ask user h-value, Width, & Height
% dlg_title = 'Potential Distribution';
% num_lines = 1;
% val = inputdlg(prompt,dlg_title,num_lines); % assign h-value

h =0.5 %str2double(val{1}); % Converting to number
diwid =10 %str2double(val{2});
diheight =20 %str2double(val{3});
dimenheight = diheight/h; % calculating number of rows for rectangle geometry
dimenwidth = diwid/h; % calculating number of columns for rectangle geometry
i = dimenheight-1; % calculating number of nodes
j = dimenwidth-1; % calculating number of nodes
nodes = i*j; % number of nodes
air = 2.997925*10^8; % Velocity of Light in Air
miu0 = 4*pi*10^(-7); % Permeability Free-Space
ep0 = 8.854*10^(-12); % Dielectric Constant (Permittivity)
m = 1; % Which Mode for m
n = 1; % Which Mode for n
rhs = 0; 
lhs = 0;
top = 0;
bot = 0;
anafreq = (air/(2*pi))*(sqrt(((m*pi/diwid)^2)+((n*pi/diheight)^2)));

% Creating A Matrix

for m = 1:j-1
    for n = 1:i
        d1a2(m+(n-1)*j, m+(n-1)*j+1) = 1; % Creating Diagonal of Ones Above -4 Matrix
        d1a2(m+(n-1)*j+1, m+(n-1)*j) = 1; % Creating Diagonal of Ones Below -4 Matrix
    end
end
for m = 1:j
    for n = 1:i-1
        d3a4(m+(n-1)*j, m+n*j) = 1; % Creating Diagonal of Ones Above Far
        d3a4(m+n*j, m+(n-1)*j) = 1; % Creating Diagonal of Ones Below Far
    end
end
A = -4*eye(nodes,nodes)+d1a2+d3a4; % Creating Matrix of Coefficients

[v c] = eig(A); % Calculate Eigenvalues & vectors of A
eigval = diag(c); % Extracting Eigenvalues
cutfreq = (1/(2*pi))*(air/h)*sqrt(-eigval); % Convert to Potential Cut off Frequencies
dif = abs(anafreq-cutfreq) % Finding which is closest to Analytical Frequencies
loc = find(abs(dif-min(dif))<0.001) % Find where lowest Difference is
cutfreq1 = cutfreq(loc(1,:),1) % Extract Closest to Ananlytical CutFreq
beta = (2*pi*anafreq)*(1/air)*sqrt(1-(cutfreq1/anafreq)) % Calc Beta 
vel = air/sqrt(1-(cutfreq1/anafreq)^2) % Calc Velocity of wave
coeff = ((2*pi*cutfreq1)^2)*miu0*ep0;

A1 = ((coeff*h^2)-4)*eye(nodes,nodes)+d1a2+d3a4;



% Calculate wave impedance

eta = 120*pi*sqrt(1-(cutfreq1/anafreq)^2)

phi = v(:,loc(1,:))'; % Calculating Phi Values
NEW = reshape(phi,j,i); % Recreates Matrix of Phis to fix Dimensions of i x j
MAP = zeros(i+2,j+1); % Empty Matrix of Zeros to Hold Next Values
for k = 1:i
    MAP(k+1, j+2) = rhs; % Creates Given Potentials in Image
    MAP(k+1, 1) = lhs;
    MAP(1,1:j+2) = top;
    MAP(i+2,1:j+2) = bot; %if there's bot value
end
x = 0:i-1; % Determines how many Rows are Needed
y = 0:j-1; % Determines of many Columns are Needed
MAP(2+x, 2+y) = NEW.'; % Matrix of Phis
%fprintf('cutoff frequency=%d',cutfreq);
fprintf('Ananlytical CutFreq=%d\n',anafreq);
fprintf('Closest Freq to Ananlytical CutFreq=%d\n',cutfreq1);

imagesc(NEW.')
contourf(NEW.'/eta)

colorbar