%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Computation for Transmission Line Geometry    %%
%% Anthony Gasbarro                              %%
%% EE 675                                        %%    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
clearvars;
close all;
clc;
%% Given Data Parameters %%

%Dimensions of TL
    height = 2; %overall height given as 2.02, rounded down to 2
    width= 7;  %overall width

%Conductor Parameters
    widthConductor = 1.00; 
    thicknessConductor = 0; %Assume thickness to be zero
%Dielectric Parameters
    heightDielectric = 1.00; %dielectric height per given data
    ErDielectric = 9.6;
    
%Universal Constants
    E0 = 8.854*10^(-12); %permittivity of free space (F/m)
    c = 299792458; %speed of light in a vaccuum(m/s)


%% Adjustable parameters %%

%potential of the conductor 
    phiConductor = 10.00; 

%size of mesh
    h = 0.1;
    
%Contour Distance from Center Conductor
%this is a value that the distance from the center conductor is divided by
    contour = 2;
    
%User can designate potential for edges
    leftPotential = 0;
    rightPotential = 0;
    topPotential = 0;
    bottomPotential = 0;

%% Generation of Matrices A and B %%

%dimensions of points for phi matrix
    colMax = (width/h)-1;
    rowMax = (height/h)-1; 

%determine 'n' for A matrix    
    nodesTotal = rowMax*colMax;
    
%determine given parameter dimensions as nodes in the phi matrix
    
    %height of dielectric nodes along y axis, and section above
        yDielectric = (heightDielectric/h);
        yAboveDielectric = rowMax - yDielectric;
        
    %width of conductor nodes along x axis, and sections to side
        xConductorWidth = (widthConductor/h) + 1;
        xNextToConductor = (colMax - xConductorWidth)/2;
    

%The whole computation will run two times consecutively
%the first run will find capacitance with all air 
%the second run will compute capacitance with Dielectric
for run = 1:2   
    if(run == 1)
        Er = 1; %set capacitance for air
    else
        Er = ErDielectric;
        capacitanceAir = capacitance;
        outputMatrixAir = outputMatrix;
        phiMatrixAir = phiMatrix;
        AmatrixAir = Amatrix;
        BmatrixAir = Bmatrix;
    end
    
    % generate 'A' Matrix with diagonal -4's, along with empty 'B' matrix
    Amatrix = eye(nodesTotal, nodesTotal)*-4;
    Bmatrix = zeros(nodesTotal, 1);
   
    % iterate through each row to fill in 'A' and 'B' matrices
    % 'i' is used as the iterator throughout the 'A' and 'B' matrices
    % x,y are used to determine position within the phi matrix
    for i = 1:nodesTotal

        %create an x,y point for reference in the phi matrix
        %y is determined as the modulus of x since, y-axis for tables 
        %are inversed from typical i,j notation
        x = mod(i,colMax);
        y = fix(i/colMax)+1; 

        %account for rows that are divisible by rmax
        if(x == 0)
            x = colMax;
            y = y-1;
        end

        %Conditional flags for special conditions
        %above, below, left, right of, and at conductor
        %left, right, top, bottom edges of mesh
        %conditional flags are set by if statements below 
        %all flags are initialized to zero at each loop iteration

            aboveConductor = 0; %current x-position is just above conductor
            belowConductor = 0; %current x-position is just below conductor
            onConductor = 0; %current x,y position is on top of the conductor
            leftOfConductor = 0;
            rightOfConductor = 0;
            leftEdge = 0;
            rightEdge = 0;
            topEdge = 0;
            bottomEdge = 0;

        % The following computations are broken into two major sections
        % when at the dielectric interface and when not at the interface
        % Base case is when not at dielectric interface
        if(y ~= yAboveDielectric+1)

            %Determine if at nodes above/below conductor
            %first check if x-position is within conductor bounds
            %then check if at a node directly above or below and mark flag
            %then finally calculate value stored in Bmatrix by subtracting
            if((x > xNextToConductor) && (x <= (xNextToConductor + xConductorWidth)))
                if(y == yAboveDielectric)
                    aboveConductor = 1;
                    Bmatrix(i) = Bmatrix(i) - phiConductor;
                elseif(y == yAboveDielectric + 2)
                    belowConductor = 1;
                    Bmatrix(i) = Bmatrix(i) - phiConductor;
                end
            end 

            % This section determines if at an Edge in the phi matrix
            % and subtracts potential from 'B' matrix element 
            % first a flag is set to determine if at a boundary,
            % then potential is updated based on input values above  
            if(x == 1)
                leftEdge = 1;
                Bmatrix(i) = Bmatrix(i) - leftPotential;
            end

            if(x == colMax)
                rightEdge=1;
                Bmatrix(i) = Bmatrix(i) - rightPotential;
            end

            if(y == 1)
                topEdge=1;
                Bmatrix(i) = Bmatrix(i) - topPotential;
            end

            if(y == rowMax)
                bottomEdge=1;
                Bmatrix(i) = Bmatrix(i) - bottomPotential;
            end

            %The set the 1's if no edge or conductor nearby
            if(leftEdge == 0)
                Amatrix(i,i-1) = 1;
            end
            if(rightEdge == 0)
                Amatrix(i,i+1) = 1;
            end
            if(topEdge == 0 && belowConductor == 0)
                Amatrix(i, i-colMax) = 1;
            end
            if(bottomEdge == 0 && aboveConductor == 0)
                Amatrix(i, i+colMax) = 1;
            end


        % The 'else' case below handles second major section where 
        % 'y' is AT THE DIELECTRIC INTERFACE
        % this section is broken into the following cases: 
        %   -between conductor and edge 
        %   -on the conductor 
        %   -directly to the left/right of conductor,
        %   -and at an matrix edge 
        %   -if at interface with no edge or conductor
        else

            % Determine if NOT on the conductor and compute potential
            if(x <= xNextToConductor || x>xConductorWidth+xNextToConductor)
                Amatrix(i,i) = -4*(Er + 1);

            % if ON conductor set flag and set 'A' and 'B' matrix values  
            else
                onConductor = 1;
                Amatrix(i,i) = 1;
                Bmatrix(i) = phiConductor;

            end

            % The following section accounts for if directly to the left then if
            % directly to the right of the conductor 
            if(onConductor == 0)
                if(x == xNextToConductor) 
                    leftOfConductor = 1;
                    Bmatrix(i) = Bmatrix(i) - (Er + 1)*phiConductor;

                elseif(x == (xConductorWidth + xNextToConductor + 1))
                    rightOfConductor = 1;
                    Bmatrix(i) = Bmatrix(i) - (Er + 1)*phiConductor;  
                end

                % The next section is for determining edge calculations 
                % while at the dielectric interface
                if(x == 1)
                    leftEdge = 1;
                    Bmatrix(i) = Bmatrix(i) - (Er + 1)*leftPotential;
                end

                if(x == colMax)
                    rightEdge = 1;
                    Bmatrix(i) = Bmatrix(i) - (Er + 1)*rightPotential;
                end

                if(y == 1)
                    topEdge = 1;
                    Bmatrix(i) = Bmatrix(i) - topPotential*2;
                end 

                if(y == rowMax)
                    bottomEdge = 1;
                    Bmatrix(i) = (i) - Er*bottomPotential*2;
                end

                %This section handles cases where there is no conductor and
                %no edge adjacent to the node
                if(leftEdge == 0 && rightOfConductor == 0) 
                    Amatrix(i, i-1) = Er+1; 
                end
                if(rightEdge == 0 && leftOfConductor == 0)
                    Amatrix(i, i+1) = Er+1;
                end
                if(topEdge == 0)
                    Amatrix(i, i-colMax) = 2;
                end
                if(bottomEdge == 0)
                    Amatrix(i, i+colMax) = Er*2;
                end
            end
        end
    end

    %% Generate phi Matrix and place values in Output matrix


    % Multiply 'A' matrix by inverse 'B' matrix to create Phi Matrix
    phiMatrix = mldivide(Amatrix, Bmatrix);

    % Create Output Matrix to store Phi Outputs
    outputMatrix = zeros(rowMax+2, colMax+2);

    %Create dimension references for output matrix
        rowMaxOutput = rowMax + 2;
        colMaxOutput = colMax + 2;

        % set user input potential values to output matrix boundaries
        for i = 2:rowMax+1
            outputMatrix(i,1) = leftPotential;
        end

        for i = 2:rowMax+1
            outputMatrix(i,rowMax+2) = rightPotential;
        end

        for i = 2:colMax+1
            outputMatrix(1, i) = topPotential;
        end

        for i = 2:colMax+1
            outputMatrix(rowMax+2, i) = bottomPotential;
        end

        % Iterate through phi matrix and place values inside of output Matrix
        i = 1;
        for x = 2:rowMax+1
            for y = 2:colMax+1
                 outputMatrix(x,y) = phiMatrix(i);
                 i = i+1;
            end
        end

     %% Compute Capacitance Based on Output Matrix and Contour   

    % initialize charge value to zero
        q = 0;

    % The following calculates the contour based on the fraction input when
    % initalizing variables
        leftContour = fix((colMaxOutput - xConductorWidth)/(2*contour))+1;
        rightContour = colMaxOutput - leftContour;
        topContour = fix((rowMaxOutput - yDielectric)/contour)+1;
        bottomContour = rowMaxOutput - fix(yDielectric/contour);

     yDielectricInterface =  rowMaxOutput-yDielectric - 1; 
     
     %Iterate through contour and determine charge by adding up values in each
     %cell of the output matrix three situations are taken into account whether
     %a particular node is at the above, at, or below the dielectric interface
     %first iterate top to bottom doing both inside and outside left/right sides
     for i = topContour+1:bottomContour-1

         %for above the Dielectric Interface
         %only factor E0, permitivity in air
         %take both sides and divide by two
         if(i < yDielectricInterface)
             q = q +(E0*(outputMatrix(i, leftContour-1)- outputMatrix(i, leftContour+1))/2);
             q = q +(E0*(outputMatrix(i, rightContour+1)-outputMatrix(i, rightContour-1))/2);

         %for when AT the Dielectric Interface
         %factor both E0 and ER and divide by 1/4th
         elseif(i == yDielectricInterface)
             q = q +(E0*(outputMatrix(i, leftContour-1)- outputMatrix(i, leftContour+1))/4);
             q = q +(E0*(outputMatrix(i, rightContour+1)-outputMatrix(i, rightContour-1))/4);
             q = q +(Er*E0*(outputMatrix(i, leftContour-1)- outputMatrix(i, leftContour+1))/4);
             q = q +(Er*E0*(outputMatrix(i, rightContour+1)-outputMatrix(i, rightContour-1))/4);
         else
             q = q +(Er*E0*(outputMatrix(i, leftContour-1)- outputMatrix(i, leftContour+1))/2);
             q = q +(Er*E0*(outputMatrix(i, rightContour+1)-outputMatrix(i, rightContour-1))/2);

         end

     end

    %next iterate from left to right while computing top and bottom
    %top side only factors E0 since it is above Dielectric Interface
    for i = leftContour+1:rightContour-1
        q = q +(E0*(outputMatrix(topContour-1, i)- outputMatrix(topContour+1, i))/2);
        q = q +(Er*E0*(outputMatrix(bottomContour+1, i)-outputMatrix(bottomContour-1, i))/2);
    end

    %use final accumulated charge value to compute capacitance
    capacitance = -q/phiConductor;

end
%% Compute Impedance and Velocity of Propegation

Z = 1/(c*sqrt(capacitanceAir*capacitance));
Vp = c*sqrt(capacitanceAir/capacitance);

%% Generate output and Display on Figures

%Display Air output on heatmap
figure('Name','Air Run','NumberTitle','off');
imagesc(outputMatrixAir)
hold on
colormap(hot)
colorbar
hold off

%Display Dielectric output on heatmap
figure('Name','Dielectric Run','NumberTitle','off');
imagesc(outputMatrix)
hold on
colormap(hot)
colorbar
hold off

fprintf("*************Computation Complete*************\n")
fprintf("Capacitance with Air = %d\n",capacitanceAir);
fprintf("Capacitance with Dielectric = %d\n",capacitance);
fprintf("Characteristic Impedance = %d\n", Z);
fprintf("Velocity of Propegation = %d\n",Vp);