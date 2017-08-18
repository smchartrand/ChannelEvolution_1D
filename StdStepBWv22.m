function [ yout,counter ] = StdStepBWv22(  B,yd,S,Sfd,j,Qw,Dgx,dx,Nn,counter,HLimd,HLimu  )
    
    %% STANDARD STEP METHOD: CODE TO COMPUTE WATER SURFACE PROFILE
    % This code solves the gradually varied flow equation using the standard 
    % step method. The calculations begin by assuming an initial water depth.
    % This particular function is set-up for backwater calculations.
    
    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION
    % Gravitational acceleration - meters per square second
    g = 9.81; 
    % Two times g
    gg = 2 * g;
    
    % Tolerance for iterative solution
    tol = 0.001;
    
    %% SETUP THE CALCULATION TO ONLY ENTER THE WHILE LOOP FOR NODES 
    % DOWNSTREAM OF THE MOST UPSTREAM NODE.
    if counter == 1 || yd(1,j+1) == 0
        
        if j == Nn
            
            % Assign the boundary condition depth at node j+1 as the limiting
            % depth for Fr > 1 (meters)
            ybc = HLimd(1,j);
            Bbc = B(1,Nn);
            
        else
            
            ybc = HLimu(1,j);
            Bbc = B(1,j+1);
            
        end
        
        % Compute area (square meters)
        Abc = B(1,j) * ybc;
        % Compute wetted perimeter (meters)
        Pbc = B(1,j) + (2 * ybc);
        % Trial Mannings n
        Manningsbc = 0.042 * ((Dgx(1,j) * 3.281)^(1/6));
        % Trial k parameter for the friction slope calculation
        kbc = (Manningsbc * Qw) ^ 2;
        % Compute the predicted friction slope
        Sfbc = kbc * Pbc ^ 1.333 / Abc ^ 3.333;
        % Set the initial value of dh
        dh = 1;
        
        % Pass the ycrit value to a different parameter to avoid issues
        ypass = ybc + ((S(1,j) - Sfbc) * (-dx));

        %% BEGIN STANDARD STEP CALCULATION

        while abs(dh) > tol

            % Step 1: Assign a trial depth to start the calucations.  Trial depth 
            % is equivalent to the normal depth.   
            ytrial = ypass;

            % Step 2: Compute the trial area, wetted perimeter, etc. based on the 
            % trial depth.  Trial Mannings is based on the Strickler relation.
            % Trial area
            Atrial = B(1,j) * ytrial;
            % Trial wetted perimeter
            Ptrial = B(1,j) + (2 * ytrial);
            % Trial Mannings n
            Mannings = 0.042 * ((Dgx(1,j) * 3.281)^(1/6));
            % Trial k parameter for the friction slope calculation
            ktrial = (Mannings(1,j) * Qw) ^ 2;
           
            % Step 3: Compute the trial head based on steps 1 and 2.
            headtrial = ytrial + ((Qw / Atrial)^2 / gg);

            % Step 4: Compute the friction slope for the present node based on the
            % estimates from steps 1 - 3.
            Sftrial = ktrial * Ptrial ^ 1.333 / Atrial ^ 3.333;

            % Step 5: Compute the friction head loss based on results from steps 1
            % - 4 and from subcritical results or supercritical results if
            % at counter >= 2
            headlosstrial = (S(1,j) - ((Sfbc + Sftrial) / 2)) * (-dx);

            % Step 6: Compute second trial head based on results of steps 1 - 5.
            Areac = ybc * Bbc;
            headtrial2 = (ybc + ((Qw / Areac) ^ 2 / gg)) + headlosstrial;

            % Step 7: Compute the difference between the assumed (headtrial) and 
            % the computed (headtrial2).  When the difference is less than the 
            % threshold break the loop. 
            if headtrial2 - headtrial > tol
                
                % Increment the looping parameter to avoid infinite loop
                ypass = ypass + tol;
                
            elseif headtrial - headtrial2 > tol
                
                ypass = ypass - tol;
            
            end
            
            dh = headtrial2 - headtrial;
            
            if ypass < HLimd(1,j)
                
                ypass = HLimd(1,j);
                % break the loop if at the limiting Froude depth
                break
                
            end
            
        end
        
    elseif counter > 1
    
        dh = 1;
        
        % Pass the ycrit value to a different parameter to avoid issues
        ypass = yd(1,j+1) + ((S(1,j) - Sfd(1,j+1)) * (-dx));

        %% BEGIN STANDARD STEP CALCULATION

        while abs(dh) > tol

            % Step 1: Assign a trial depth to start the calucations.  Trial depth 
            % is equivalent to the normal depth.   
            ytrial = ypass;

            % Step 2: Compute the trial area, wetted perimeter, etc. based on the 
            % trial depth.  Trial Mannings is based on the Strickler relation.
            % Trial area
            Atrial = B(1,j) * ytrial;
            % Trial wetted perimeter
            Ptrial = B(1,j) + (2 * ytrial);
            % Trial Mannings n
            Mannings = 0.042 * ((Dgx(1,j) * 3.281)^(1/6));
            % Trial k parameter for the friction slope calculation
            ktrial = (Mannings * Qw) ^ 2;
           
            % Step 3: Compute the trial head based on steps 1 and 2.
            headtrial = ytrial + ((Qw / Atrial)^2 / gg);

            % Step 4: Compute the friction slope for the present node based on the
            % estimates from steps 1 - 3.
            Sftrial = ktrial * Ptrial ^ 1.333 / Atrial ^ 3.333;

            % Step 5: Compute the friction head loss based on results from steps 1
            % - 4 and from subcritical results or supercritical results if
            % at counter >= 2
            headlosstrial = (S(1,j) - ((Sfd(1,j+1) + Sftrial) / 2)) * (-dx);

            % Step 6: Compute second trial head based on results of steps 1 - 5.
            Area = (yd(1,j+1) * B(1,j+1));
            headtrial2 = (yd(1,j+1) + ((Qw / Area) ^ 2 / gg)) + headlosstrial;

            % Step 7: Compute the difference between the assumed (headtrial) and 
            % the computed (headtrial2).  When the difference is less than the 
            % threshold break the loop. 
            if headtrial2 - headtrial > tol
                
                % Increment the looping parameter to avoid infinite loop
                ypass = ypass + tol;
                
            elseif headtrial - headtrial2 > tol
                
                ypass = ypass - tol;
            
            end
            
            dh = headtrial2 - headtrial;
            
            if ypass < HLimd(1,j)
                
                ypass = HLimd(1,j);
                % break the loop if at the limiting Froude depth
                break
                
            end
            
        end

    end
    
    % Pass the final value for the converged depth to the yu parameter
    yout(1,j) = ypass;
    dhout(1,j) = dh;
        
end

