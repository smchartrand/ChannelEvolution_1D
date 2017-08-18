function [ yout ] = PredCorrTestFDG( j,yd,S,Scrit,dx,qw,Ksx )

    %% MIXED FLOW APPROXIMATION: CODE TO COMPUTE RVF WATER SURFACE PROFILE
    % This code solves the backwater equation according to a predictor -
    % corrector scheme using a limiting Fr condition coupled with a 
    % limiting depth and solves the forewater equation according to the 
    % stnadrad step numerical scheme. The backwater numerical scheme is 
    % sourced from Parker's e-book Chapter 20 and the forewater numerical
    % scheme is sourced from Henderson 1966 and Chow 1963.

    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION
    % Gravitational acceleration - meters per square second
    g = 9.81; 
    % Manning-Strickler alpha parameter - dimensionaless
    a = 8.1;
    
    %% STEP 2: RUN THE LOOP BACKWARD TO COMPUTE THE BACKWATER PROFILE  
                              
        % PREDICTOR STEP
        % Compute the predicted incremental depth change
        Fbackpd(1,j) = (((S(1,j+1) / Scrit(1,j+1)) - FrSqd(1,j+1)) / (1 - FrSqd(1,j+1))) * Scrit(1,j+1);
        % Compute predicted water depth at current location
        ypd(1,j) = yd(1,j+1) - Fbackpd(1,j) * dx;
        % Compute the square of the predicted Froude number
        FrSqpd(1,j) = qw(1,j)^2 / g / (ypd(1,j)^3);
        Frpd(1,j) = sqrt(FrSqpd(1,j));
        % Compute the predicted Chezy friction coefficient
        Cfpd(1,j) = (1/a^2) * (ypd(1,j)/Ksx(1,j))^(-1/3);
        % Compute the predicted friction slope
        Sfpd(1,j) = Cfpd(1,j) * FrSqpd(1,j);

        % CORRECTOR STEP
        % Compute corrected incremental depth change
        Fbackcd(1,j) = (((S(1,j) / Scrit(1,j)) - FrSqd(1,j)) / (1 - FrSqpd(1,j))) * Scrit(1,j);
        % Compute average water depth slope
        dyd(1,j) = (Fbackpd(1,j) + Fbackcd(1,j)) / 2;
        % Compute water depth based on the limiting condition
        yout(1,j) = yd(1,j+1) - (dyd(1,j)*dx);    
            
end



