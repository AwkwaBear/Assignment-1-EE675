% EE 471 | Ex. 4.17

% Khaldoon Ismaeel 
% 2017 
%
% ------------------------------------------------------------------------
%
% This program is designed to use the Method of Moments to find the
% capacitance as a function of the seperation of distance between parallel
% plates.
%
%
% 
% ------------------------------------------------------------------------
clear
clc


% Variables
W = 1; % width of capacitor plate
d = 1; % distance between plates
Area = W*W; % area of plate
n = 15;
N = n*n; % # of sections
l = W/n; % length of subarea
SA = l*l; % subarea
pi = 3.14159;           %value of pi
DA = pi*(l/2)^2;
V = 1; % potential difference between plates
ep0 = 8.854*10^(-12); % value of epsilon
inc = .1;
distinc = d/inc;
diag = sqrt(DA)/(2*ep0*sqrt(pi)); % diag values  (2*sqrt(pi)
coef = SA/(4*pi*ep0); % coef values for other subareas of plate
rtop = zeros(n);
rbot = zeros(n);
seta = zeros(2*N,2*N);
count = 1;
%counter = 1;
stop = n*n;
% cap = zeros(1,N);
q = zeros(1,distinc);
disvec = zeros(2*N);
var = 1;
counting = 1;
% Creating Matrix of Distances
for counter = 1:distinc
    for varx = 1:n
        for vary = 1:n
            for i = 1:n
                for j = 1:n
                    rtop(i,j) = sqrt(((abs(vary-i))*l)^2 + ((abs(varx-j))*l)^2);
                    rbot(i,j) = sqrt((counting*inc)^2 + (rtop(i,j))^2);
                end
            end
        Rtop = rtop(:)';
        Rbot = rbot(:)';
        Row1 = horzcat(Rtop,Rbot);
        Row2 = horzcat(Rbot,Rtop);
        if count <= N
            disvec(count,:) = Row1;
            disvec(N+count,:) = Row2;
            count = count+1;
        end
        end
    end
    if counting < distinc
        count = 1;
    end
    vec = disvec + eye(2*N);
    B = 1./vec; % Distance is in denominator
    A2 = (ones(2*N)*coef) - eye(2*N)*coef; % Calculating the Coef for not diag
    A3 = eye(2*N)*diag; % Coeff for diag
    A4 = A2 + A3; % Add two togehter
    tote = A4.*B; % Complete Matrix
    bpos = ones(N,1); % vector of ones
    bneg = -ones(N,1); % vector of -ones
    btot = vertcat(bpos,bneg);
    rho = tote\btot; % Calculating total rho Values from +
    halfpos = rho(1:N);
    halfneg = rho(N+1:2*N); % Calculating total rho Values from -
    R1 = reshape(halfpos,n,n);
    R0 = -R1;
    q(1,counting) = sum(halfpos*SA);
    counting = counting+1;
end

dist = inc:inc:d;
dwratio = dist./W;
prod = 1./dwratio.*W;
C0 = prod.*W.*ep0; % analytical capacitance
cap = q./(V-(-V));

capratio = cap./C0;

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
