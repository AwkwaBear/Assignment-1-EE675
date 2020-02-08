%EE471  Transmission Line 
clear;
clc;

%user-defined variables:
width = 3;          %width of cross section in m
height = 2;         %height of cross section in m
h = 0.1;            %mesh size in m
dieheight = 1;		%bottom dielectric's height in m
conwidth = 1;		%conductor strip width in m
                    %(note: conwidth/h should be even/odd if width/h is even/odd)
                    %(note: assume conductor height is negligible)
diepsr = 9.6;       %dielectric's relative permittivity
                    %(note: top medium is air)
conphi = 10;        %conductor's potential in V
eps0 = 8.854*10^(-12); %permittivity of free space in F/m
v0 = 299792458;     %speed of light in vacuum in m/s
l = 0;              %left edge potential in V
r = 0;              %right edge potential in V
t = 0;              %top edge potential in V
b = 0;              %bottom edge potential in V
airphi = 0;         %potential in air
contour = 2;        %set contour size, must be greater than 1


%x by y unknown points inside the rectangle
x = (width/h) - 1;
y = (height/h) - 1;
%dielectric height in terms of y
diey = dieheight/h;
%conductor width in terms of x
conx = (conwidth/h) + 1;
%points along the y axis above the dielectric
abovediey = y - diey;
%points along x axis on each side of the conductor
sideconx = (x - conx)/2;

%The below if statement checks that all parameters are legal before running. 
if (h>0 && mod(width,h)==0 && mod(width,h)==0 && mod(dieheight,h)==0 && mod(conwidth,h)==0 && mod(width-conwidth,h*2)==0 && width>h && height>h && dieheight>h && conwidth>=h && diey<=y-1 && conx<=x-2)

    %will run twice, first to find C with all air, second for C with dielectric
    for caprun = 1:2
        if(caprun==1)
            epsr=1;
        elseif(caprun==2)
            epsr=diepsr;
            cap0=cap;
            R0=R;
            phi0=phi;
            A0=A;
            B0=B;
        end
        
        %create matrix A and vector B for central difference method equation
        %in A*phi = B form
        A = zeros(x*y,x*y);
        B = zeros(x*y,1);
        
        %iterate through cross section's inner points
        %starting from left to right, top to bottom,
        for i = 1:x*y
            
            %find position along x axis and inverted y axis
            xposition = mod(i,x);
            yposition = fix(i/x)+1;
            if(xposition==0)
                xposition = x;
                yposition = yposition-1;
            end
            
            %conditional flags
            edgeleft=0;
            edgeright=0;
            edgetop=0;
            edgebottom=0;
            conleft=0;
            conright=0;
            contop=0;
            conbottom=0;
            atcon=0;
                      
            %if not at dielectric interface
            if(yposition~=abovediey+1)
                %subtract 4 times potential at center point
                A(i,i) = -4;
                
                %if directly above or below conductor
                if(xposition>sideconx && xposition<=sideconx+conx)
                    %if conductor next on bottom
                    if(yposition==abovediey)
                        B(i) = B(i) - conphi;
                        conbottom=1;
                        %if conductor next on top
                    elseif(yposition==abovediey+2)
                        B(i) = B(i) - conphi;
                        contop=1;
                    end
                end
                
                %subtract edge potential normally from Laplace equation if next to edge
                %if next to left edge
                if(xposition==1)
                    B(i) = B(i) - l;
                    edgeleft=1;
                end
                %if next to right edge
                if(xposition==x)
                    B(i) = B(i) - r;
                    edgeright=1;
                end
                %if next to top edge
                if(yposition==1)
                    B(i) = B(i) - t;
                    edgetop=1;
                end
                %if next to bottom edge
                if(yposition==y)
                    B(i) = B(i) - b;
                    edgebottom=1;
                end
                
                %add unknown neighbor potential normally if there is no edge or conductor
                %if no left edge
                if(edgeleft==0)
                    A(i,i-1) = 1;
                end
                %if no right edge
                if(edgeright==0)
                    A(i,i+1) = 1;
                end
                %if no top edge and no top conductor
                if(edgetop==0 && contop==0)
                    A(i,i-x) = 1;
                end
                %if no bottom edge and no bottom conductor
                if(edgebottom==0 && conbottom==0)
                    A(i,i+x) = 1;
                end
                
            else %if at dielectric interface
                %if between conductor and left/right edge
                if(xposition<=sideconx || xposition>sideconx+conx)
                    %subtract 4 times (epsilon_r + 1) times potential at center point
                    A(i,i) = -4*(epsr+1);
                    
                else %otherwise at conductor
                    A(i,i) = 1;
                    B(i) = conphi;
                    atcon=1;
                end
                
                %if not at conductor, perform operations
                if(atcon==0)
                    %if conductor next on right
                    if(xposition==sideconx)
                        B(i) = B(i) - (epsr+1)*conphi;
                        conright=1;
                        %if conductor next on left
                    elseif(xposition==sideconx+conx+1)
                        B(i) = B(i) - (epsr+1)*conphi;
                        conleft=1;
                    end
                    
                    %if also next to edges while at dielectric interface
                    %use interface equation
                    %if next to left edge
                    if(xposition==1)
                        B(i) = B(i) - (epsr+1)*l;
                        edgeleft=1;
                    end
                    %if next to right edge
                    if(xposition==x)
                        B(i) = B(i) - (epsr+1)*r;
                        edgeright=1;
                    end
                    %if next to top edge
                    if(yposition==1)
                        B(i) = B(i) - 2*t;
                        edgetop=1;
                    end
                    %if next to bottom edge
                    if(yposition==y)
                        B(i) = B(i) - 2*epsr*b;
                        edgebottom=1;
                    end
                    
                    %add unknown neighbor potential if there is no edge or conductor
                    %use interface equation
                    %if no left edge and no left conductor
                    if(edgeleft==0 && conleft==0)
                        A(i,i-1) = epsr+1;
                    end
                    %if no right edge and no right conductor
                    if(edgeright==0 && conright==0)
                        A(i,i+1) = epsr+1;
                    end
                    %if no top edge
                    if(edgetop==0)
                        A(i,i-x) = 2;
                    end
                    %if no bottom edge
                    if(edgebottom==0)
                        A(i,i+x) = 2*epsr;
                    end
                end
                
            end
            
        end
        
        
        %solve system of equations
        phi = A\B;
        
        
        %organize potentials graphically on rectangle R
        R = zeros(y+2,x+2);
        
        %insert arbitrary corner values
        corner = 0;
        R(1,1) = corner;
        R(1,x+2) = corner;
        R(y+2,1) = corner;
        R(y+2,x+2) = corner;
        
        %left
        for n = 2:y+1
            R(n,1)=l;
        end
        %right
        for n = 2:y+1
            R(n,x+2)=r;
        end
        %top
        for n = 2:x+1
            R(1,n)=t;
        end
        %bottom
        for n = 2:x+1
            R(y+2,n)=b;
        end
        
        n=1;
        for i = 2:y+1
            for j = 2:x+1
                R(i,j)= phi(n);
                n=n+1;
            end
        end
        
        %find capacitance using phi values in R matrix
        %initialize charge value
        q=0;
        %get x by y dimensions of R matrix
        newy=y+2;
        newx=x+2;
        %set capacitance contour edges about halfway between conductor and edge
        topcapy = fix((newy-diey)/contour)+1;
        bottomcapy = newy-fix(diey/contour);
        leftcapx = fix(((newx-conx)/2)/contour)+1;
        rightcapx = newx-leftcapx;
        
        %find charge based on contour
        %going from top to bottom along right and left sides of contour
        for i = topcapy+1:bottomcapy-1
            %if above dielectric interface
            if(i<newy-diey-1)
                q=q+eps0*(R(i,leftcapx-1)-R(i,leftcapx+1))/2;
                q=q+eps0*(R(i,rightcapx+1)-R(i,rightcapx-1))/2;
                %if at dielectric interface
            elseif(i==newy-diey-1)
                q=q+eps0*(R(i,leftcapx-1)-R(i,leftcapx+1))/4;
                q=q+epsr*eps0*(R(i,leftcapx-1)-R(i,leftcapx+1))/4;
                q=q+eps0*(R(i,rightcapx+1)-R(i,rightcapx-1))/4;
                q=q+epsr*eps0*(R(i,rightcapx+1)-R(i,rightcapx-1))/4;
                %if below dielectric interface
            else
                q=q+epsr*eps0*(R(i,leftcapx-1)-R(i,leftcapx+1))/2;
                q=q+epsr*eps0*(R(i,rightcapx+1)-R(i,rightcapx-1))/2;
            end
        end
        %going from left to right along top and bottom sides of contour
        for i = leftcapx+1:rightcapx-1
            %for top side
            q=q+eps0*(R(topcapy-1,i)-R(topcapy+1,i))/2;
            %for bottom side
            q=q+epsr*eps0*(R(bottomcapy+1,i)-R(bottomcapy-1,i))/2;
        end
        
        %calculate capacitance
        cap = -q/conphi;
        
    end
    
    %find impedance Z0
    Z0 = 1/(v0*sqrt(cap*cap0));
    %find velocity of propagation Vp
    Vp = v0*sqrt(cap0/cap);
    
    
    %display various values (warning: small mesh size may take long to display)
%     
%     A0
%     A
%     B0
%     B
%     phi0
%     phi
%     R0
%     R
    
    %create blank spaces to show contour on plot
    for i = topcapy:bottomcapy
        %for left side
        R(i,leftcapx)=0;
        %for right side
        R(i,rightcapx)=0;
    end
    for i = leftcapx:rightcapx
        %for top side
        R(topcapy,i)=0;
        %for bottom side
        R(bottomcapy,i)=0;
    end
    
    %display graphical rectangle
    figure(1)
    imagesc(R)
    hold on
    colormap(jet)
    colorbar
    figure(2) 
    imagesc(R0)
    colormap(jet)
    colorbar
%     quiver (x,y)
    %display parameters
    fprintf('Successfully Finished, Chosen Parameters: \n')
    fprintf('%.2f m by %.2f m cross section with %.4f m mesh size \n', width, height, h)
    fprintf('%.2f m high dielectric with %.2f m wide conductor centered on it \n', dieheight, conwidth)
    fprintf('Relative Permittivity of Dielectric: %.2f, Conductor Potential: %.2f V \n', epsr, conphi)
    fprintf('Edge Potentials: Left: %.2f V, Right: %.2f V, Top: %.2f V, Bottom: %.2f V \n\n', l, r, t, b)
    fprintf('Resulting Parameters: \n')
    fprintf('Capacitance with Dielectric: %d, Capacitance with All Air: %d \n', cap, cap0)
    fprintf('Characteristic Impedance: %3.2f, Velocity of Propagation: %d \n', Z0, Vp)

else
    fprintf('Error - Inappropriate Dimensions: \n')
    fprintf('%.2f m by %.2f m cross section with %.4f m mesh size \n', width, height, h)
    fprintf('%.2f m high dielectric with %.2f m wide conductor centered on it \n', dieheight, conwidth)
end