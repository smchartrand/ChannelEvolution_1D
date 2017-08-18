function [ yu,WSEu,Tbedeu,Vu,Ushearu,Fru,Sfu ] = HydraulicJump( B,dx,Ksx,Dgx,Qw,qw,S,hj,yd,yu,WSEu,Tbedeu,Vu,Ushearu,Fru,Sfu,n1 )

    %% UNTITLED Summary of this function goes here
    % Detailed explanation goes here
    
    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION
    % Relative density - dimensionless
    Rd = (2650/1000) - 1;             
    % Water density - kg per cubic meter
    Rho = 1000;
    % Gravitational acceleration - meters per square second
    g = 9.81; 
    % Manning-Strickler alpha parameter - dimensionaless
    a = 8.1;
    
    hjindex = find(hj ~= 0);
    L = length(hjindex);
    hjindexfirst = min(hjindex);
    hjindexlast = max(hjindex);

    for i = L:-1:1

        %% STANDARD STEP METHOD: CODE TO COMPUTE WATER SURFACE PROFILE
        % This code solves the gradually varied flow equation using the standard 
        % step method. The calculations begin by assuming an initial water depth.
        % This particular function is set-up for backwater calculations.

        %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION
        % Tolerance for iterative solution
        tol = 0.001;

        %% PRE-ALLOCATE LOOPING PARAMETERS AND OTHER PARAMETERS
        ytrial(1,L)=0;
        Atrial(1,L)=0;
        Ptrial(1,L)=0;
        Rtrial(1,L)=0;
        Mannings(1,L)=0;
        ktrial(1,L)=0;
        Sftrial(1,L)=0;
        headtrial(1,L)=0;
        headtrial2(1,L)=0;
        headlosstrial(1,L)=0;
        yuu(1,L)=0;

        %% SETUP THE CALCULATION TO ONLY ENTER THE WHILE LOOP FOR NODES 
        % DOWNSTREAM OF THE MOST UPSTREAM NODE.
        if i == L

            % Assign the boundary condition depth at node j+1 as the limiting
            % depth for Fr > 1 (meters)
            ybc = yd(1,(hjindex(i)+1));
            Bbc = B(1,hjindex(i));

            % Compute area (square meters)
            Abc = B(1,hjindex(i)) * ybc;
            % Compute wetted perimeter (meters)
            Pbc = B(1,hjindex(i)) + (2 * ybc);
            % Trial Mannings n
            Manningsbc = 0.042 * ((Dgx(1,hjindex(i)) * 3.281)^(1/6));
            % Trial k parameter for the friction slope calculation
            kbc = (Manningsbc * Qw) ^ 2;
            % Compute the predicted friction slope
            Sfbc = kbc * Pbc ^ 1.333 / Abc ^ 3.333;

            dh = 1;

            % Pass the ycrit value to a different parameter to avoid issues
            ypass = ybc + ((S(1,hjindex(i)) - Sfbc) * (-dx));

            %% BEGIN STANDARD STEP CALCULATION

            while abs(dh) > tol

                % Step 1: Assign a trial depth to start the calucations.  Trial depth 
                % is equivalent to the normal depth.   
                ytrial(1,i) = ypass;

                % Step 2: Compute the trial area, wetted perimeter, etc. based on the 
                % trial depth.  Trial Mannings is based on the Strickler relation.
                % Trial area
                Atrial(1,i) = B(1,hjindex(i)) * ytrial(1,i);
                % Trial wetted perimeter
                Ptrial(1,i) = B(1,hjindex(i)) + (2 * ytrial(1,i));
                % Trial hydraulic radius
                Rtrial(1,i) = Atrial(1,i) / Ptrial(1,i);
                % Trial Mannings n
                Mannings(1,i) = 0.042 * ((Dgx(1,hjindex(i)) * 3.281)^(1/6));
                % Trial k parameter for the friction slope calculation
                ktrial(1,i) = (Mannings(1,i) * Qw) ^ 2;

                % Step 3: Compute the trial head based on steps 1 and 2.
                headtrial(1,i) = ytrial(1,i) + ((Qw / Atrial(1,i))^2 / (2 * g));

                % Step 4: Compute the friction slope for the present node based on the
                % estimates from steps 1 - 3.
                Sftrial(1,i) = ktrial(1,i) * Ptrial(1,i) ^ 1.333 / Atrial(1,i) ^ 3.333;

                % Step 5: Compute the friction head loss based on results from steps 1
                % - 4 and from subcritical results or supercritical results if
                % at counter >= 2
                headlosstrial(1,i) = (S(1,hjindex(i)) - ((Sfbc + Sftrial(1,i)) / 2)) * (-dx);

                % Step 6: Compute second trial head based on results of steps 1 - 5.
                headtrial2(1,i) = (ybc + ((Qw / (ybc * Bbc)) ^ 2 / (2 * g))) + headlosstrial(1,i);

                % Step 7: Compute the difference between the assumed (headtrial) and 
                % the computed (headtrial2).  When the difference is less than the 
                % threshold break the loop. 
                if headtrial2(1,i) - headtrial(1,i) > tol

                    % Increment the looping parameter to avoid infinite loop
                    ypass = ypass + tol;

                elseif headtrial(1,i) - headtrial2(1,i) > tol

                    ypass = ypass - tol;

                end

                dh = headtrial2(1,i) - headtrial(1,i);

            end

        elseif i < L

            dh = 1;

            % Pass the ycrit value to a different parameter to avoid issues
            ypass = yuu(1,i+1) + ((S(1,hjindex(i)) - Sfuu(1,i+1)) * (-dx));

            %% BEGIN STANDARD STEP CALCULATION

            while abs(dh) > tol

                % Step 1: Assign a trial depth to start the calucations.  Trial depth 
                % is equivalent to the normal depth.   
                ytrial(1,i) = ypass;

                % Step 2: Compute the trial area, wetted perimeter, etc. based on the 
                % trial depth.  Trial Mannings is based on the Strickler relation.
                % Trial area
                Atrial(1,i) = B(1,hjindex(i)) * ytrial(1,i);
                % Trial wetted perimeter
                Ptrial(1,i) = B(1,hjindex(i)) + (2 * ytrial(1,i));
                % Trial hydraulic radius
                Rtrial(1,i) = Atrial(1,i) / Ptrial(1,i);
                % Trial Mannings n
                Mannings(1,i) = 0.042 * ((Dgx(1,hjindex(i)) * 3.281)^(1/6));
                % Trial k parameter for the friction slope calculation
                ktrial(1,i) = (Mannings(1,i) * Qw) ^ 2;

                % Step 3: Compute the trial head based on steps 1 and 2.
                headtrial(1,i) = ytrial(1,i) + ((Qw / Atrial(1,i))^2 / (2 * g));

                % Step 4: Compute the friction slope for the present node based on the
                % estimates from steps 1 - 3.
                Sftrial(1,i) = ktrial(1,i) * Ptrial(1,i) ^ 1.333 / Atrial(1,i) ^ 3.333;

                % Step 5: Compute the friction head loss based on results from steps 1
                % - 4 and from subcritical results or supercritical results if
                % at counter >= 2
                headlosstrial(1,i) = (S(1,hjindex(i)) - ((Sfuu(1,i+1) + Sftrial(1,i)) / 2)) * (-dx);

                % Step 6: Compute second trial head based on results of steps 1 - 5.
                headtrial2(1,i) = (yuu(1,i+1) + ((Qw / (yd(1,i+1) * B(1,i+1))) ^ 2 / (2 * g))) + headlosstrial(1,i);

                % Step 7: Compute the difference between the assumed (headtrial) and 
                % the computed (headtrial2).  When the difference is less than the 
                % threshold break the loop. 
                if headtrial2(1,i) - headtrial(1,i) > tol

                    % Increment the looping parameter to avoid infinite loop
                    ypass = ypass + tol;

                elseif headtrial(1,i) - headtrial2(1,i) > tol

                    ypass = ypass - tol;

                end

                dh = headtrial2(1,i) - headtrial(1,i);

            end
            
        end
        % Pass the final value for the converged depth to the yu parameter
        yuu(1,i) = ypass;
        % Compute water surface elevation
        WSEuu(1,i) = n1(1,hjindex(i)) + yu(1,i);
        % Trial area
        Auu(1,i) = B(1,hjindex(i)) * yuu(1,i);
        % Trial wetted perimeter
        Puu(1,i) = B(1,hjindex(i)) + (2 * yuu(1,i));
        % Trial Mannings n
        Mannings(1,i) = 0.042 * ((Dgx(1,hjindex(i)) * 3.281)^(1/6));
        % Trial k parameter for the friction slope calculation
        kuu(1,i) = (Mannings(1,i) * Qw) ^ 2;
        % Compute the predicted friction slope
        Sfuu(1,i) = kuu(1,i) * Puu(1,i) ^ 1.333 / Auu(1,i) ^ 3.333;
        % Compute the velocity based on the contiuity equation
        Vuu(1,i) = qw(1,hjindex(i)) / (yuu(1,i));
        % Compute dimensional shear stress distribution per node to node slope
        Tbeduu(1,i) = Rho * Vuu(1,i)^2 * (a^-2 * (yuu(1,i) / Ksx(1,hjindex(i)))^-0.33);
        % Compute the shear velocity
        Ushearuu(1,i) = sqrt(Tbeduu(1,i) / Rho);
        % Compute dimensionless shear stress distribution per node to node slope
        Tbedeuu(1,i) = Ushearuu(1,i)^2 / (g * Rd * Dgx(1,hjindex(i)));
        % Compute the square of the predicted Froude number
        FrSquu(1,i) = qw(1,hjindex(i))^2 / g / (yuu(1,i)^3);
        Fruu(1,i) = sqrt(FrSquu(1,i));

    end
    % Swap out water depth values just computed for the hydraulic
    % jump region for those of the yu vector using logical indexing
    yu([hjindexfirst hjindexlast]) = yuu(:);
    WSEu([hjindexfirst hjindexlast]) = WSEuu(:);
    Tbedeu([hjindexfirst hjindexlast]) = Tbedeuu(:);
    Vu([hjindexfirst hjindexlast]) = Vuu(:);
    Ushearu([hjindexfirst hjindexlast]) = Ushearuu(:);
    Fru([hjindexfirst hjindexlast]) = Fruu(:);
    Sfu([hjindexfirst hjindexlast]) = Sfuu(:);    

end

