%EE675 Assignment Method of Moments Antenna Problem by Khaldoon Ishmael
%Reference on pg767 Example 9.4 

clear all
clc

% Variables:
c = 2.998*(10^8);        % speed of light                    
f = 1 * 10^9;            % Frequency in GHz
eta0 = 120*pi;           %impedance in free space eta0=377
eps0 = 8.854*10^(-12);   %permittivity of free space in F/m
mu0 = 4*pi*10^(-7);      % permeability of free space in N/(A^2)
lambda = (2*pi*c)/f;     % wavelength lambda in m

% Case 1: L1=0.94, l/lambda=1/2
% Case 2: L1=0.47, l/lambda=1/4
% Case 3: L1=0.705, l/lambda=3/8
% Case 4: L1=0.23, l/lambda<<1

L1= 0.94;              % Antenna length in m
a1= 0.0001;             % Antenna radius in m 

L = lambda*L1;
a = lambda*a1;
% a=L/518;

n=60;                   % total number of segments

Vg=1;                    % voltage difference between gap
pi = 3.14159;           %value of pi
B0 = (2*pi)/lambda;   % phase constant beta0


dz = L/n;                                % length of segments
z = (-L/2 + dz/2):dz:(L/2 - dz/2);       % antenna centered at 0, finds center of each segment
z(n+1) = L/2 - dz;

% choose arbitrary additional point to vector of z
    % this point would be placed just under the top point
  
    % Create B matrix 
B = ((1j*(Vg/(2*eta0)))*sin(B0*abs(z)))';

%create matrix A and vector B in A*J = B form
 A = zeros(n+1,n+1);
    
    
    for o = 1:n+1   %"point o" is observation point
                    % "point : " is source point
        %+1 to solve for C at point slightly away from center of a section
        %z = 0; %m at center of antenna,z = n*dl; %m at one rod end
    R = sqrt(a^2 + (z-z(o)).^2);  
                      % z is the testing point, z(o) is the z'
                      
    A(o,:) = (exp(-1*1j*B0*R)./R)*dz;
                      % using the delta expansion function, we
                      % can eliminate the integral.  This is applied to
                      % each row
    end
       
A(:,n+1) = cos(B0*z)';  % last column of A matrix

current = (A)\B*4*pi;
current = current(1:n);
antenna_length = z(1:n);

%Calculation
current_center = current(abs(n/2));
Impedance = Vg/(current_center)

%display plot
plot(abs(current),antenna_length/lambda)
title('Current Density Distribution on Antenna')
ylabel('Section of Antenna-Both Rods(m)') % y-axis label
xlabel('Current Density (A/(m^2))') % x-axis label
 
%display Results.
fprintf('Successfully Finished :-) \n')
fprintf('Length: %.5f m \nRadius: %.5f m \nSections: %.8f m\n', L, a, dz)
fprintf('Gap Potential: %.2f V \nWavelength: %.2f m\n\n', Vg, lambda)

