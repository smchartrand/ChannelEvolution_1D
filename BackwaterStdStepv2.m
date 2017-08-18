function [ Tbede,Ushear,y,WSE,V,S ] = BackwaterStdStepv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn )

    %% MIXED FLOW APPROXIMATION: CODE TO COMPUTE RVF WATER SURFACE PROFILE
    % This code solves the backwater equation according to a predictor -
    % corrector scheme using a limiting Fr condition coupled with a 
    % limiting depth and solves the forewater equation according to the 
    % stnadrad step numerical scheme. The backwater numerical scheme is 
    % sourced from Parker's e-book Chapter 20 and the forewater numerical
    % scheme is sourced from Henderson 1966 and Chow 1963.

    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION
    % Relative density - dimensionless
    Rd = (2650/1000) - 1;             
    % Water density - kg per cubic meter
    Rho = 1000;
    % Gravitational acceleration - meters per square second
    g = 9.81; 
    % Manning-Strickler alpha parameter - dimensionaless
    a = 8.1;
    % Momentum coefficient - similar to alpha - dimensionless
    Beta = 1.0;
    
    %% PRE-ALLOCATE MEMORY FOR PARAMETERS
    S(1,Nn)=0;
    Sfd(1,Nn)=0;
    FrSqd(1,Nn)=0;
    Frd(1,Nn)=0;
    yd(1,Nn)=0;
    ynorm(1,Nn)=0;
    WSEd(1,Nn)=0;
    Vd(1,Nn)=0;
    Ad(1,Nn)=0;
    Pd(1,Nn)=0;
    HRd(1,Nn)=0;
    kd(1,Nn)=0;
    Taud(1,Nn)=0;
    Tstarg(1,Nn)=0;
    Ustard(1,Nn)=0;
    SpFod(1,Nn)=0;
    Mannings(1,Nn)=0;
    Head(1,Nn)=0;
    nc(1,Nn)=0;
    ncindex(1,Nn)=0;
    yout(1,Nn)=0;
    yr(1,Nn)=0;
    order(1,Nn)=0;
    counter2=0;
    HC(1,Nn)=0;
    
    %% FIRST COMPUTE SOME INITIAL TERMS TO GET THINGS STARTED.
    % Compute the critical water depth (m)
    ycrit = (qw.^2 ./ (g)).^(1/3);
    % Compute the critical bed slope
    Scrit = (Ksx.^0.33 .* qw.^2) ./ (a^2 * g * (ycrit.^(1/0.30)));
    % Set the limiting Froude condition for backwater calculations
    FrLimd = 0.9999999999999999999;
    % Compute the limiting flow depth associated with the limiting lower Froude
    HLimd = (qw ./ g^0.5 ./ FrLimd) .^(2/3);
    % Set the limiting Froude condition for backwater calculations
    FrLimu = 1.0000000000000000001;
    % Compute the limiting flow depth associated with the limiting upper Froude
    HLimu = (qw ./ g^0.5 ./ FrLimu) .^(2/3);
        
    %% STEP 1: LOOP THROUGH AND COMPUTE THE NODE TO NODE BED SLOPE AND
    % the normal depth. Also flag locations that transition from a pool to a riffle, 
    % or a deep to a shallow. This equates to finding locations where the derivative 
    % field changes direction.
    for j = 1:1:Nn
       
        % If loop is at the upstream most node
        if j == 1
            % Compute the bed slope at the u/s boundary
            S(1,j) = (n1(1,j) - n1(1,j+1)) / dx;
        
        % If loop is at all stations in between the boundaries
        elseif j > 1 && j < Nn
            % Compute the bed slope at all other stations
            S(1,j) = ((n1(1,j) - n1(1,j+1)) / dx);
            
        % If loop is at all other stations
        else
            % Compute the bed slope at the d/s boundary
            S(1,j) = (n1(1,j-1) - n1(1,j)) / dx;
            
        end
        
        % Compute the normal depth approximation if the bed slope is > 0
        if S(1,j) > 0 
            ynorm(1,j) = ((Ksx(1,j)^0.33*qw(1,j)^2)/(a^2*g*S(1,j)))^0.30;
        
        % Else set the normal depth to the appropriate limiting depth
        else
            ynorm(1,j) = HLimd(1,j);
            
        end
        
        % Call a function to compute the sign of the bed elevation
        % derivative for use in flagging pool locations
        [ order,nc,ncindex ] = DerivativeSignv2( j,n1,Nn,dx,order,nc,ncindex );
        
    end
    
    % Specify the pool slope from which to compute the residual depth. The 
    % pool slope is set as the mean bed slope at all stations.
    Smp = mean(S);
    
    % Compute a hydraulic index to track whether conditions are mild or
    % steep.  HC > 1 indicate mild conditions and HC < 1 indicate steep.
    HC = abs(Scrit ./ S);
    
    % Call the Derivative function to flag locations where the
    % derivative field changes directions, and then assign a control 
    % if it applies. This could be done in a different loop but I will do 
    % it here for now.
    
    %ncindex = find(nc~=0);
    ncsum = sum(ncindex);
    
    if ncsum ~= 0

        [ yr ] = ControlElevationv2( n1,Nn,nc,ncindex,yr,S );
        
    end

    % Clear the looping parameter for use next
    clear j
    
    %% STEP 2: RUN THE LOOP BACKWARD TO COMPUTE THE BACKWATER PROFILE
    % Initialize a counter used in the standard step calculations
    
    for j = Nn:-1:1
        % Direct the profile calcuations to the first node in the
        % subcritical profile where the computed depth is less than the
        % critical depth    
            
            if j == Nn 
                % Initialize the first counter
                counter = 1;
                
            else
                % Advance the counter
                counter = counter + 1;

            end
            
            if S(1,j) > 0
                
                if j == Nn
                    
                    if HC(1,j) > 1
                        
                        yout(1,j) = ynorm(1,j);
                        %yout(1,j) = DWSE - n1(1,end);
                        
                    else
                        
                        yout(1,j) = ynorm(1,j);
                        %yout(1,j) = DWSE - n1(1,end);
                    
                    end
                    
                else
                    
                    [ yout,counter ] = StdStepBWv22( B,yd,S,Sfd,j,Qw,Dgx,dx,Nn,counter,HLimd,HLimu );
                    
                    if yout(1,j) < HLimd(1,j)

                        yout(1,j) = HLimd(1,j);

                    end
                        
                end
                
            else
                
                [ yout,counter2] = AdvSlopeRCM( yout,j,counter2,yr,dx,HLimd,Smp );
                
            end
            
            % Assigning the standard step solution for depth to the yu vector
            yd(1,j) = yout(1,j);              
            % Compute water surface elevation
            WSEd(1,j) = n1(1,j) + yd(1,j);
            if j < Nn
                if WSEd(1,j) < WSEd(1,j+1)
                    WSEd(1,j) = WSEd(1,j+1) + (dx * tan(Smp));
                    yd(1,j) = WSEd(1,j) - n1(1,j);
                end
            end
            % Compute area (square meters)
            Ad(1,j) = B(1,j) * yd(1,j);
            % Compute wetted perimeter (meters)
            Pd(1,j) = B(1,j) + (2 * yd(1,j));
            % Compute hydraulic radius (meters)
            HRd(1,j) = Ad(1,j) / Pd(1,j);
            % Trial Mannings n
            Mannings(1,j) = 0.042 * ((Dgx(1,j) * 3.281)^(1/6));
            % Trial k parameter for the friction slope calculation
            kd(1,j) = (Mannings(1,j) * Qw) ^ 2;
            % Compute the velocity based on the contiuity equation
            Vd(1,j) = qw(1,j) / (yd(1,j));
            % Compute the total head for the computed flow depth
            Head(1,j) = yd(1,j) + (((Qw / Ad(1,j))^2 / (2 * g)));
            % Compute the shear velocity - meters per second
            Ustard(1,j) = ((Ksx(1,j)^(1/3) * (qw(1,j)^2)) / (a^2))^(3/20) * g^(7/20) * S(1,j)^(7/20);
            % Compute friction coefficient from Manning-Strickler formulation
            Cfd(1,j) = (1/a^2) * (yd(1,j)/Ksx(1,j))^(-1/3);
            % Compute the dimensional stress based on gradually varied flow - Pascals
            Taud(1,j) = (Cfd(1,j) * Rho * Vd(1,j)^2);
            % Compute the dimensionaless shear stress Shields number
            % according to the Shields equation. I tried the
            % Manning-Strickler formulation but it gave very large
            % estimates of dimensionaless stress
            Tstarg(1,j) = Taud(1,j) / (Rho * Rd * g * Dgx(1,j));
            % Compute the square of the predicted Froude number
            FrSqd(1,j) = qw(1,j)^2 / g / (yd(1,j)^3);
            Frd(1,j) = sqrt(FrSqd(1,j));
            % Compute the predicted friction slope
            Sfd(1,j) = kd(1,j) * Pd(1,j) ^ 1.333 / Ad(1,j) ^ 3.333;
            % Compute the specific force (momentum per Henderson, 1966)
            SpFod(1,j) = ((Qw * Beta) / (g * B(1,j) * yd(1,j))) + ((B(1,j) * yd(1,j)) * (yd(1,j) / 2));
 
    end

    clear j    
         
    %% STEP 4: USE RESULTS OF THE BACKWATER CALCULATIONS TO 
    % determine the appropriate profile solution. Backwater solution is
    % used at nodes where the backwater specific momentum governs;
    % forewater solution is used at nodes where the forewater specific
    % momentum governs. Try to handle this using logical indexing to
    % optimize for speed.
    % Write the computed profile results to the proper out parameters
    y = yd;
    WSE = WSEd;
    Tbede = Tstarg;
    Tbed = Taud;
    V = Vd;        
    Ushear = Ustard;
    %Fr = Frd; 
    %Cf = Cfd;
            
end



