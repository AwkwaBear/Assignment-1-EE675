% EE471, Khaldoon Ismaeel/ Rod

clc;
clear;
user1 = 'How many division do you want in the rod? \n';
section  = input(user1);

%user-defined variables:

v=1;                    %potential at center of section in V
l = 1;                    %length of cylindrical rod in m
r = 0.001;               %radius of cylindrical rod in m 
e = 8.854*10^(-12);     %permittivity of free space in F/m
pi = 3.14159;            %value of pi
mu0 = 4*pi*10^(-7);     %permeability of free space in N/(A^2)

% basic equations
h = 1/section;  %calculate h sections of rod based on section width
nodx = l/h+1;  %length within section for finding diagonal element

if (h>0)&& (l>h)  %check if dimensions appropriate

% Matrix A
A = zeros(nodx);
for i = 1:nodx  %iterate through rod's sections
    for j = 1:nodx 
        A(i,j) = (2*pi*r*h)/abs(i-j);
        
    end
end
for i = 1:nodx
    A(i,i) = (4*pi*r*log(h/r));
end

% Matrix B
B = ones(nodx,1)*4*pi*e; %set voltage at center of section and multiply by 4*pi*e from formula

% inverse calculations of the two matrixes
C = A\B;

%display various values (warning: small mesh size may take long to display)
     %  A
    %    B
    %    C

charge = C*2*pi*r*l;

% Axis value
Y = [charge;charge(numel(charge))];
X = (0:1/(numel(charge)):1);

% Plot of the data
%plot(X,Y)
stairs(X,Y);
title('Charges of a Finite Rod');
xlabel('Length [m]');
ylabel('charges [1/ro]');

%display parameters
fprintf('Successfully Finished, Chosen Parameters: \n')
fprintf('Cylindrical Rod with %.5f m length and %.5f m radius divided into %10f m sections\n', l, r, h)
fprintf('Center Potential: %.2f V\n\n', v)
else  
fprintf('Error - Inappropriate Dimensions: \n')
end 
