function [ Tbede,Ushear,y,WSE,S ] = MixedFlowTestv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn )

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
    Tbedd(1,Nn)=0;
    Tbeded(1,Nn)=0;
    Usheard(1,Nn)=0;
    SpFod(1,Nn)=0;
    
    Sfu(1,Nn)=0;
    FrSqu(1,Nn)=0;
    Fru(1,Nn)=0;
    yu(1,Nn)=0;  
    WSEu(1,Nn)=0;
    Vu(1,Nn)=0;
    Au(1,Nn)=0;
    Pu(1,Nn)=0;
    HRu(1,Nn)=0;
    Mannings(1,Nn)=0;
    ku(1,Nn)=0;
    Tbedu(1,Nn)=0;
    Tbedeu(1,Nn)=0;
    Ushearu(1,Nn)=0;
    SpFou(1,Nn)=0;
    Head(1,Nn)=0;
    nc(1,Nn)=0;
    ncindex(1,Nn)=0;
    hj(1,Nn)=0;
    yout(1,Nn)=0;
    yr(1,Nn)=0;
    order(1,Nn)=0;
    counter2=0;
    
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
    % Specify a minimum pool slope from which to compute the residual depth
    Smp = 0.001;
        
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
    
    % Call the Derivative function to flag locations where the
    % derivative field changes directions, and then assign a control 
    % if it applies. This could be done in a different loop but I will do 
    % it here for now.
    
    %ncindex = find(nc~=0);
    ncsum = sum(ncindex);
    
    if ncsum ~= 0

        [ nc,yr ] = ControlElevationv2( n1,Nn,nc,ncindex,yr,S );
        
    end

    % Clear the looping parameter for use next
    clear j
    
    %% STEP 2: RUN THE LOOP BACKWARD TO COMPUTE THE BACKWATER PROFILE
    % Initialize a counter used in the standard step calculations
    
    for j = Nn:-1:1
        % Direct the profile calcuations to the first node in the
        % subcritical profile where the computed depth is less than the
        % critical depth
        if ynorm(1,j) < ycrit(1,j)           
            % Step into the backwater calculation function. The backwater
            % calculation is set-up according to either the standard step
            % method for cases of S > 0 and the flow-depth gradient GVF
            % equation for cases of S <= 0.  The later has been coded as a
            % predictor-corrector scheme.
            if j == Nn || ynorm(1,j+1) > ycrit(1,j+1)
                % Initialize the first counter
                counter = 1;
                
            else
                % Advance the counter
                counter = counter + 1;

            end
            
            if S(1,j) > 0
                
                [ yout,counter ] = StdStepBWv2( B,yd,S,Sfd,j,Qw,Dgx,dx,Nn,counter,HLimd,HLimu );
                
            else
                
                [ yout,counter2] = AdvSlopeRCM( yout,j,counter2,yr,dx,HLimu,Smp );
                
            end
            
            % Assigning the standard step solution for depth to the yu vector
            yd(1,j) = yout(1,j);              
            % Compute water surface elevation
            WSEd(1,j) = n1(1,j) + yd(1,j);
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
            % Compute dimensional shear stress distribution per node to node slope
            Tbedd(1,j) = Rho * Vd(1,j)^2 * (a^-2 * (yd(1,j) / Ksx(1,j))^(-1/3));
            % Compute the shear velocity
            Usheard(1,j) = sqrt(Tbedd(1,j) / Rho);
            % Compute dimensionless shear stress distribution per node to node slope
            Tbeded(1,j) = Usheard(1,j)^2 / (g * Rd * Dgx(1,j));
            % Compute the square of the predicted Froude number
            FrSqd(1,j) = qw(1,j)^2 / g / (yd(1,j)^3);
            Frd(1,j) = sqrt(FrSqd(1,j));
            % Compute the predicted friction slope
            Sfd(1,j) = kd(1,j) * Pd(1,j) ^ 1.333 / Ad(1,j) ^ 3.333;
            % Compute the specific force (momentum per Henderson, 1966)
            SpFod(1,j) = ((Qw * Beta) / (g * B(1,j) * yd(1,j))) + ((B(1,j) * yd(1,j)) * (yd(1,j) / 2));

        end
 
    end

    clear j    
         
    %% STEP 3: RUN THE LOOP FORWARD TO COMPUTE THE FOREWATER PROFILE.
    % First use logical indexing to flag locations where the slope is
    % greater than the critical slope
    
    for j = 1:1:Nn
        
        % Direct the profile calcuations to the first node in the
        % subcritical profile where the computed depth is less than the
        % critical depth
        if ynorm(1,j) > ycrit(1,j)
            % Step into the forewater calculation function. The forewater
            % calculation is set-up according to the standard step method.
            if j == 1 || ynorm(1,j-1) < ycrit(1,j-1)
                % Initialize the first counter
                counter = 1;

            else
                % Advance the counters
                counter = counter + 1;
                
            end
            
            % Step into the forewater calculation function. The forewater
            % calculation is set-up according to the standard step method.
            % The uses the Newton-Raphson method to solve the non-uniform
            % flow equation.
            [ yout,counter ] = StdStepFWv2( B,yu,S,Sfu,j,Qw,Dgx,dx,Nn,counter,HLimu );
            
            % Assigning the standard step solution for depth to the yu vector
            yu(1,j) = yout(1,j);          
            % Compute water surface elevation
            WSEu(1,j) = n1(1,j) + yu(1,j);
            % Compute area (square meters)
            Au(1,j) = B(1,j) * yu(1,j);
            % Compute wetted perimeter (meters)
            Pu(1,j) = B(1,j) + (2 * yu(1,j));
            % Compute hydraulic radius (meters)
            HRu(1,j) = Au(1,j) / Pu(1,j);
            % Trial Mannings n
            Mannings(1,j) = 0.042 * ((Dgx(1,j) * 3.281)^(1/6));
            % Trial k parameter for the friction slope calculation
            ku(1,j) = (Mannings(1,j) * Qw) ^ 2;
            % Compute the velocity based on the contiuity equation
            Vu(1,j) = qw(1,j) / (yu(1,j));
            % Compute the energy head for the computed flow depth
            Head(1,j) = yu(1,j) + (((Qw / Au(1,j))^2 / (2 * g)));
            % Compute dimensional shear stress distribution per node to node slope
            Tbedu(1,j) = Rho * Vu(1,j)^2 * (a^-2 * (yu(1,j) / Ksx(1,j))^-0.33);
            % Compute the shear velocity
            Ushearu(1,j) = sqrt(Tbedu(1,j) / Rho);
            % Compute dimensionless shear stress distribution per node to node slope
            Tbedeu(1,j) = Ushearu(1,j)^2 / (g * Rd * Dgx(1,j));
            % Compute the square of the predicted Froude number
            FrSqu(1,j) = qw(1,j)^2 / g / (yu(1,j)^3);
            Fru(1,j) = sqrt(FrSqu(1,j));
            % Compute the predicted friction slope
            Sfu(1,j) = ku(1,j) * Pu(1,j) ^ 1.333 / Au(1,j) ^ 3.333;
            % Compute the specific force (momentum per Henderson, 1966)
            SpFou(1,j) = ((Qw * Beta) / (g * B(1,j) * yu(1,j))) + ((B(1,j) * yu(1,j)) * (yu(1,j) / 2));
            % Flag instances when the flow depth is getting greater. This
            % signifies the occurrence of a hydraulic jump. When the end of
            % a hydraulically steep sub-reach is encountered, step into the
            % hydraulic jump function to estimate the flow depth through
            % the hydraulic jump. Re-write all the pertinent values in the
            % relevant vectors.
            if counter > 1 && yu(1,j) > yu(1,j-1)
                
                hj(1,j) = 1;
                
            end  
           
        end
        
    end
    
    clear j
    
    %% STEP 4: USE RESULTS OF FOREWATER FLAGGING FOR LOCATIONS WHERE FLOW
    % depth increased toward critical depth (flagged above with the hj
    % parameter) to recompute the profile for a hydraulic jump or
    % subcritical regime.
    
%     if yu(1,j+1) == 0
%                    
%         [ yu,WSEu,Tbedeu,Vu,Ushearu,Fru,Sfu ] = HydraulicJump( B,dx,Ksx,Dgx,Qw,qw,S,hj,yd,yu,WSEu,Tbedeu,Vu,Ushearu,Fru,Sfu,n1 );
% 
%     end
       
    %% STEP 4: USE RESULTS OF THE BACKWATER AND FOREWATER CALCULATIONS TO 
    % determine the appropriate profile solution. Backwater solution is
    % used at nodes where the backwater specific momentum governs;
    % forewater solution is used at nodes where the forewater specific
    % momentum governs. Try to handle this using logical indexing to
    % optimize for speed.
    % Write the computed profile results to the proper out parameters
    yu(yu == 0) = yd(yd > 0);
    y = yu;
    WSEu(WSEu == 0) = WSEd(WSEd > 0);
    WSE = WSEu;
    Tbedeu(Tbedeu == 0) = Tbeded(Tbeded > 0);
    Tbede = Tbedeu;
    Vu(Vu == 0) = Vd(Vd > 0);
    V = Vu;        
    Ushearu(Ushearu == 0) = Usheard(Usheard > 0);
    Ushear = Ushearu;
    Fru(Fru == 0) = Frd(Frd > 0);
    Fr = Fru; 
            
end



