function [  ] = BedEvolutionFunction( varargin )

    %% BEDEVOLUTION: CODE TO COMPUTE BED EVOLUTION ALONG 1-DIMENSIONAL CHANNEL
    % This code computes bed evolution based on user specified hydraulic 
    % and sediment transport schemes.  Code can proceed with normal flow 
    % or backwater flow approximation, uniform or mixed grain size sediment 
    % mixtures, and simple excess stress or hiding function based sediment
    % transport formulations. Suspended sediment is ignored and the model 
    % has been explicitely constructed to explore simplifications of real
    % streams. Much of the model is informed by Parker's e-book Chapters 4,
    % 7, 14, 17 and 20.

    format long e
    clc
    close all
    
    % Time to completion status bar 
    hh = waitbar(1,'Model Initialization');
    
    % Find the waitbar object
    hw=findobj(hh,'Type','Patch');
    
    % Changes the color to green
    set(hw,'EdgeColor',[0 0 1],'FaceColor',[0 0 1]) % changes the color to green

    %% LOAD AN INPUT DATA FILE IN -ASCII FORMAT TO SET-UP THE MODEL
    % The file can be printed to PDF for documentation purposes of model
    % conditions.
    [Input1] = load('BedEvolutionInput1_MH.txt');
    [Input2] = load('BedEvolutionInput2_MH.txt');
    [Input3] = load('BedEvolutionInput3_MH.txt');
    [Input3b] = load('BedEvolutionInput3b_MH.txt');
    [Input4] = load('BedEvolutionInput4_MH.txt');
    
    %% DEFINE HYDRAULIC, SEDIMENT TRANSPORT, HYDROLOGIC AND MOBILITY CONDITIONS
    % Define hydraulic approximation: '1' for Normal Flow, '2' for
    % Backwater or '3' for Mixed Flow (not working yet).
    HYD = Input1(18,1);
    
    % Define sediment transport approximation: '1' for Excess Stress, 
    % '2' for Wilcock-Crowe, '3' for Ashida-Michue, or '4' for Ashida-Michue applied with the 
    % friction angle mobility model.
    SED = Input1(19,1);
    
    % Define if the the simulation is steady or unsteady: '1' for steady
    % and '2' for undsteady
    FLOW = Input1(20,1);
    
    % Define mobility approximation conditions: '1' for business as usual,
    % utilizing parameterizations of any sediment transport function and
    % '2' to utilize the friction angle mobility model
    MOBILITY = Input1(21,1);
    
    %% DEFINE COMPUTATIONAL SPACE, TEMPORAL AND BOUNDARY CONDITIONS
    % Spatial Conditions
    % Virtual river lenght (L) (meters)
    L = Input1(3,1);
    % Specify the station array (meters)
    Station = Input2(:,1);
    % Number of spatial nodes
    NoP = length(Station);
    % Spatial step (meters)
    dx = Input1(1,1);
    % Spatial step array (meters)
    N = Station';                
    Nn = length(N);
    % Specify the channel width (meters)
    BB = Input2(:,2);
    B = BB';
    % Specify the channel slope (fractional percent)
    SS = Input2(:,3);
    So = SS';
    % Specify and store an initial slope
    IS = SS(1,1);
    % Specify dry sediment density measured in the lab
    Density = Input1(17,1);
    % Downstream boundary condition for bed elevation calculation
    DBBC = 2;
   
    % Temporal Conditions
    % Total run time (T) (seconds)
    % If running a hydrograph make sure the total here is consistent with
    % the hydrograph data series loaded below.
    T = Input1(4,1);
    % Time step (seconds)
    dt = Input1(2,1);
    % If running steady flow build the teim step array and specify vector
    % length
    if FLOW == 1
        % Time step array (seconds)
        t = dt:dt:T;
        tt = length(t);
    end
    % Number of saves
    NoS =Input1(5,1);
    % Timestep of output saves
    timesave = (T / dt) / NoS;
    
    % Boundary Conditions
    % Water and unit water discharge - cubic and square meters per second
    Qw = Input1(6,1);
    qw = Qw ./ B;
    Qlps = Qw * 1000; %#ok<*NASGU>
    % Downstream water surface elevation boundary condition
    DWSE = Input1(16,:);
    % If running steady flow assign and convert the sediment transport rate
    if FLOW == 1
        % Unit sediment feed rate - kilograms per minute
        Qsfkgm = Input1(7,1);
        % Convert unit sediment feed rate to square meters per second
        Qsf = Qsfkgm / 60 / Density / B(1,1);
    end
    % Downstream water surface elevation
    WSbc = Input1(16,1);
    
    %% DEFINE CONSTANTS RELEVENT TO SEDIMENT TRANSPORT
    % Bed porosity - dimensionless
    Lamda = Input1(9,1);
    % Nikaradsie coefficient
    Nk = Input1(8,1);     
    % Critical dimensionless shear stress for SED = 1 - dimensionless
    Tcrit = Input1(10,1);   
    % Upwinding coefficient - from Parker: 0.5 for Normal depth, > 0.5 for Backwater
    if HYD == 1
        alpha = Input1(13,1);
    else
        alpha = Input1(14,1);
    end
    % Bed material exchange coefficient - generally taken to range from 0.2 to 0.3
    alphased = Input1(11,1);
    % Grain scale constant to compute active layer thickness - generally ranges from 1.5 to 3
    GrainScaleConstant = Input1(12,1);
    
    %% DEFINE WATER AND SEDIMENT DISCHARGE BOUNDARY CONDITIONS    
    % Downstream bed elevation to establish initial profile for Normal flow
    No = Input1(15,1);
    n11 = Input2(:,4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Refinement made for Dry Creek modeling. Check if this is a problem
    % for other model runs
    n111 = n11';
    n1 = fliplr(n111);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% PRE-ALLOCATE MEMORY FOR MODEL PARAMETERS 
    gl = length(Input3);
    Gl = gl - 1;
    pLi(1,Gl)=0;
    fi(Nn,Gl)=0;
    fix(Nn,Gl)=0;
    Psi(1,Gl)=0;
    GSm(1,Gl)=0;
    Fi(Nn,Gl)=0;
    Fix(Nn,Gl)=0;
       
    n1save1(1,Nn)=0;                         
    Ssave1(1,Nn)=0;
    Ssave(NoS,Nn)=0;
    BedSlope(NoS,Nn)=0;
    y(1,Nn)=0;
    ysave1(1,Nn)=0;
    ysave(NoS,Nn)=0;
    WSEsave1(1,Nn)=0;
    WSEsave(NoS,Nn)=0; 
    Vsave1(1,Nn)=0;
    Vsave(NoS,Nn)=0;
    Tbede(1,Nn)=0;
    DLax(1,Nn)=0;
    
    Qbsave1(1,Nn)=0;
    Qbsave(NoS,Nn)=0;
    SGrad1(NoS,Nn)=0;
    SGrad(NoS,Nn)=0;
    etasave1(NoS,Nn)=0;
    etasave(NoS,Nn)=0;
    Finew(Nn,Gl)=0;
    Fisave1(NoS,Nn,Gl)=0;
    Fisave(NoS,Nn,Gl)=0;
    
    Sigmax(1,Nn)=0;
    Sigmasave1(1,Nn)=0;
    Sigmasave(NoS,Nn)=0;
    Dgsave1(1,Nn)=0;
    Dgsave(NoS,Nn)=0;
    Dgmm(1,Nn)=0;
    DgSq1x(1,Nn)=0;
    Fssave1(1,Nn)=0;
    Fssave(NoS,Nn)=0;
    Frsave1(1,Nn)=0;
    Frsave(NoS,Nn)=0;
    AdvCounterSave(NoS)=0;
    
    Fixsum(Nn,1)=0;
    time(NoS,1)=0;
    
    D90x(1,Nn)=0;
    Dgx(1,Nn)=0;
    Dg1(1,Nn)=0;
    Lax(1,Nn)=0;
    Fsx(1,Nn)=0;
    Ksx(1,Nn)=0;
    
    %% DEFINE INITIAL GRAIN SIZE CHARACTERISTICS - 
    % BED SEDIMENT PROPERTIES FOR UNIFORM GRAIN SIZE OPTION         
    if SED == 1
        % Values are specified. Vectors need to be transposed to be
        % consistent with operations in other functions
        Dgx = Input4(:,2)';
        D90 = Input4(:,3)';
        Ksx = Nk .* D90;
    end
    
    % BED SEDIMENT PROPERTIES FOR MIXED GRAIN OPTION
    if SED == 2 || 3
        % Assign the surface grain size details to the GSD parameter
        % for use in the friction angle model.  The Input3b data is
        % transposed in order to align with the friction angle code
        GSD = Input3b';
        % Bedload texture.
        CummDistL = Input3b(end,:);    
        % Subsurface material texture.
        CummDistSu = Input3b(3:end,:);
        % Surface material texture. Three columns of data: size (mm),
        % size (psi), cummulative gradation (%)
        Grain = Input3b(1,:);
        PsiScale = Input3b(2,:);
        CummDist = Input3b(3:end,:);
        % Fill the surface grain size matrices with the correct gradations
        for j = 1:Nn
            % Compute the fractional percent of material in each grain size
            % class
            for k = 1:Gl
                % Compute the fractional percent of material in each size class bin in the upstream supplied load
                pLi(k) = CummDistL(k+1) - CummDistL(k);  
                % Compute the fractional percent of material in each size class bin in the subsurface
                fi(j,k) = CummDistSu(j,k+1) - CummDistSu(j,k); 
                % Compute the fractional percent of material in each size class bin in the surface
                Fi(j,k) = CummDist(j,k+1) - CummDist(j,k);
                % Compute characteristic grain size of each size class bin in the surface
                Psi(k) = (PsiScale(k) + PsiScale(k+1)) * 0.5;    
                % Compute the equivalent characteristic grain diameter of
                % each size class bin in the surface - meters
                GSm(k) = ((2.^ Psi(1,k))./1000);
            end
            % Assign some passing arrays to make sure the correct data is
            % used in the steps below
            CummDistpass(1,:) = CummDist(j,:);
            % Now step through the grain size statistics calculations
            fix(j,:) = fi(j,:);
            Fix(j,:) = Fi(j,:);
            CD1 = find(CummDistpass(1,:)<0.90,1,'last');
            CD2 = find(CummDistpass(1,:)>0.90,1,'first');
            D90x(j) = (2 ^ (PsiScale(1,CD1) + (0.90 - CummDistpass(1,CD1)) * ((PsiScale(1,CD2) - PsiScale(1,CD1)) / (CummDistpass(1,CD2) - CummDistpass(1,CD1))))) / 1000;
            Lax(j) = D90x(j) * GrainScaleConstant;
            Ksx(j) = Nk * D90x(j);
            Fs1 = find(Grain==2,1,'last');
            % Note the fraction sand is not multiplied by 100 - it is left
            % as a fractional number between 0 and 1.
            Fsx(j) = CummDist(j,Fs1);
            % EXECUTE THE SURFACE MATERIAL GRAINSIZESTATISTICSWORKING1 FUNCTION
            Fipass(1,:) = Fix(j,:);
            [ Sigmag,Dg ] = GrainSizeStatistics( Gl,Psi,Fipass );   
            Dgx(j) = Dg;
            Dgmm(j) = Dgx(j) * 1000;
            DgSq1x(j) = Dgmm(j)^2;
            RMSDgmm = sqrt((Nn^-1.*sum(DgSq1x)));
            Sigmax(j) = Sigmag;
        end
        clear j
    end
    
    %% EXECUTE THE HYDRUALIC MODEL - INITAL CONDITIONS
    % Data are now fed into fucntions which compute hydraulic conditions 
    % within the 1-dimensional artificial channel.  The normal function 
    % is based on uniform, steady flow.  Backwater is based on the backwater 
    % approximation for gradually varied, non-uniform flow and mixedflow is 
    % based on approximating the water surface profile through segments of 
    % supercritical and sbucritical flow.  The mixed flow function
    % presently does not work, but it is close.  The backwater function
    % can used either a predictor-corrector, standard step or RK4
    % schemes.
    if HYD == 1
        % Normal flow approximation
        [ Tbede,Ushear,y,WSE,V,S ] = Normal( n1,dx,Ksx,Dgx,Nn,qw,IS );

    elseif HYD == 2    
        % Backwater approximation to shallow water equations
        [ Tbede,Ushear,y,WSE,V,S,AdvCounter ] = BackwaterStdStepv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn );

    elseif HYD == 3    
        % Mixed-flow approximation to shallow water equations
        [ Tbede,Ushear,y,WSE ] = MixedFlowTestv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn );

    end
    
    %% EXECUTE THE FRICTION ANGLE MOBILITY MODEL - INITIAL CONDITIONS
    % Data are now fed into a function which computes probability
    % distributions of the critical Shields stress based on a froce balance
    % scaled by computed friction angles for virtual grains dropped on a
    % random bed composed of grain sizes specified by the specifed GSD of
    % the bed surface.  The model is based on the work of Wiberg and Smtih
    % (1985) and Kirchner (1990).
    if MOBILITY == 2
        % Utilized the friction angle mobility model
        [ TCMPercentileMatrix ] = FAMobilityModelSimple( GSD );
        
    end
    
    %% EXECUTE THE SEDIMENT TRANSPORT MODEL - INITIAL CONDITIONS
    % Data are now fed into functions which compute rates of sediment
    % transport.  The sed capacity funtion is based on excess shear
    % stress raised to a power and WilcockCrowe is based on percent
    % sand within the bedload material.  In the future I will add the
    % Parker and Einstein functions.  The WilcockCrowe function has
    % been optimized for speed.
    if SED == 1
        % Sediment transport capacity model
        [ Qbe ] = SedCapacity( Nn,Dgx,Tbede,Tcrit,B,Density );

    elseif SED == 2
        % Wilcock and Crowe model
        [ Qbe,PBKtranspose ] = WilcockCrowe2( Nn,Gl,Fix,Fsx,Psi,Dgx,Tbede,Ushear,Density,B );
        
    elseif SED == 3
        % Ashida-Michue model
        [ Qbe,PBKtranspose ] = AshidaMichiuev2( Nn,Gl,Fix,Dgx,Ushear,GSm );
        
    elseif SED == 4
        % Ashida-Michue model modified for a random mobility condition
        [ Qbe,PBKtranspose ] = AshidaMichiueMOB( Nn,Gl,Fix,Ushear,GSm,TCMPercentileMatrix );

    end
    
    %% SAVE THE INITIAL MODEL CONDITIONS
    % Assign an index value for data storage of initial conditions
    index = 1;

    time(index,1) = 0;
    Fisave(index,:,:) = Fi(:,:);
    Ssave(index,:) = So(1,:);
    etasave(index,:) = n1(1,:);
    Dgsave(index,:) = Dgx(1,:);
    Sigmasave(index,:) = Sigmax(1,:);
    D90save(index,:) = D90x(1,:);
    WSEsave(index,:) = WSE(1,:);
    ysave(index,:) = y(1,:);
    Vsave(index,:) = V(1,:);
    Ushearsave(index,:) = Ushear(1,:);
    if SED == 2 || 3
        RMSSize(index,:) = RMSDgmm(1,:);
    end
    Qbsave(index,:) = Qbe(1,:);    
    
    %% IF RUNNING A HYDROGRAPH NEED TO SET-UP THE HYDROGRAPH DATA SERIES
    if FLOW == 2
        % LOAD an external file specifying hydrograph for the experiment.  File
        % contains two columns of data: hydrograph time step (hours), flow rate
        % (cubic meters per second).  Hydrograph data used later in code.
        HYDRO = load('Hydrograph.txt');
        HYDTime = HYDRO(:,1) ./ dt;
        HYDFlow = HYDRO(:,2);  
        HYDLength = length(HYDTime);
        % Fill some vecotrs with zeros to manage memory
        HYDTimediff(HYDLength-1,1)=0;
        HYDFlowdiff(HYDLength-1,1)=0;
        Qwincremental(HYDLength-1,1)=0;

        % Run a loop to compute the hydrograph characteristics
        for m = 1:HYDLength-1

            HYDTimediff(m) = HYDTime(m+1) - HYDTime(m);

            HYDFlowdiff(m) =  HYDFlow(m+1) - HYDFlow(m);

            Qwincremental(m) = HYDFlowdiff(m) / HYDTimediff(m);

        end
        clear m

        % Run through the hydrograph function to generate the
        % continuous hydrograph data
        [ Qwseries ] = Hydrograph( dt,T,tt,Qw,HYDTime,Qwincremental );
        QLength = length(Qwseries);
        % Set the size of the qqw vector
        qqw(QLength,Nn)=0;

        for j = 1:QLength

            for k = 1:Nn

                qqw(j,k) = Qwseries(j,1) ./ B(k);

            end

        end
        clear j
        clear k 
    end
    
    close(hh)
    close all
    
    %% SET SOME THINGS PRIOR TO ENTERING THE MORPHODYNAMIC MODEL
    % Advance the index for data save
    index = index + 1;
        
    % Time to completion status bar 
    h = waitbar(0,sprintf('Computing - Time Simulation run time = %0.2f (hours)',T/3600));
    
    % Find the waitbar object
    hw=findobj(h,'Type','Patch');
    
    % Changes the color to green
    set(hw,'EdgeColor',[0 0 1],'FaceColor',[0 0 1]) % changes the color to green
    
    % Run time clock for the morphodynamic model
    TimeStart = tic;
    
    %% EXECUTE THE MORPHOLOGIC MODEL 
    % Temporal loop
    for i = 1:tt 
    i
        % The intial conditions are computed above so no need to repeat the
        % calculations during the first time step for hydraulics or
        % sediment transport
        if i > 1
        
            %% GRAB THE CORRECT HYDROGRAPH FLOW VALUE FOR TIMESTEP
            % This series of operations grabs the correct flow value 
            % for the given timestep and also 
            if FLOW == 2
                
                % Point to new value of flow for each time step
                Qw = Qwseries(i);
                % Unit sediment feed rate - kilograms per minute
                Qsfkgm = Qw(1,1) * 8.33333;
                % Convert unit sediment feed rate to square meters per second
                Qsf = (Qsfkgm * (1/60) * (1/Density)) / B(1,1);
                % Qsf = Qwseries(i)*8.333;
                qw = qqw(i,:);

            end          
            
            %% EXECUTE THE HYDRUALIC MODEL
            % Data are now fed into fucntions which compute hydraulic conditions 
            % within the 1-dimensional artificial channel.  The normal function 
            % is based on uniform, steady flow.  Backwater is based on the backwater 
            % approximation for gradually varied, non-uniform flow and mixedflow is 
            % based on approximating the water surface profile through segments of 
            % supercritical and sbucritical flow.  The mixed flow function
            % presently does not work, but it is close.  The backwater function
            % can used either a predictor-corrector, standard step or RK4
            % schemes.
            if HYD == 1

                [ Tbede,Ushear,y,WSE,V,S ] = Normal( n1,dx,Ksx,Dgx,Nn,qw,IS );

            elseif HYD == 2    

                [ Tbede,Ushear,y,WSE,V,S,AdvCounter ] = BackwaterStdStepv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn );

            elseif HYD == 3    

                 [ Tbede,Ushear,y,WSE ] = MixedFlowTestv2( n1,dx,B,Qw,qw,Ksx,Dgx,Nn );

            end
            
            %% EXECUTE THE FRICTION ANGLE MOBILITY MODEL - INITIAL CONDITIONS
            % Data are now fed into a function which computes probability
            % distributions of the critical Shields stress based on a froce balance
            % scaled by computed friction angles for virtual grains dropped on a
            % random bed composed of grain sizes specified by the specifed GSD of
            % the bed surface.  The model is based on the work of Wiberg and Smtih
            % (1985) and Kirchner (1990).
            if MOBILITY == 2
                % Utilized the friction angle mobility model
                [ TCMPercentileMatrix ] = FAMobilityModelSimple( GSD );

            end

            %% EXECUTE THE SEDIMENT TRANSPORT MODEL
            % Data are now fed into functions which compute rates of sediment
            % transport.  The sed capacity funtion is based on excess shear
            % stress raised to a power and WilcockCrowe is based on percent
            % sand within the bedload material.  In the future I will add the
            % Parker and Einstein functions.  The WilcockCrowe function has
            % been optimized for speed.
            if SED == 1
                % Sediment transport capacity model
                [ Qbe ] = SedCapacity( Nn,Dgx,Tbede,Tcrit,B,Density );

            elseif SED == 2
                % Wilcock and Crowe model
                [ Qbe,PBKtranspose ] = WilcockCrowe2( Nn,Gl,Fix,Fsx,Psi,Dgx,Tbede,Ushear,Density,B );
                
            elseif SED == 3
                % Ashida-Michue model
                [ Qbe,PBKtranspose ] = AshidaMichiuev2( Nn,Gl,Fix,Dgx,Ushear,GSm );
                
            elseif SED == 4
                % Ashida-Michue model modified for a random mobility condition
                [ Qbe,PBKtranspose ] = AshidaMichiueMOB( Nn,Gl,Fix,Ushear,GSm,TCMPercentileMatrix );
                
            end
       
        end
         
       %% COMPUTE THE SPATIAL SEDIMENT TRANSPORT GRADIENT AND NEW BED ELEVATIONS
       % Data are now fed into a function which computes the sediment
       % transport gradients and updated bed elevations.  The gradients are
       % based on central differences.  The function has been optimized for
       % speed based on indexing the spatial loop.
       [ sedgrad,etanew] = SedGradETANew(alpha,Qbe,Qsf,dt,dx,Lamda,n1,No,Nn,HYD,WSbc,y,DBBC);
              
        %% COMPUTE THE EXCHANGE AND SURFACE FREQUENCIES 
        if SED == 2 || SED == 3
        
            % Data are now fed into a function which computes grain
            % size distributions for each spatial node given the fractional 
            % sediment transport graidents computed in the last step.  The 
            % function has been optimized for speed by indexing the inner loop
            % and passing data into the spatial loop one node at a time.
            [ Finew,Laxp ] = ExChangeSurfaceFreqv2( sedgrad,PBKtranspose,pLi,n1,fi,Fix,Qbe,Qsf,dx,dt,Lax,DLax,i,Nn,alpha,alphased,Gl,etanew,Lamda );
            
            % COMPUTE NEW SURFACE MATERIAL STATISTICS %
            for j = 1:Nn
                
                Fixsum(j) = sum(Finew(j,:));
                l = 1:Gl;
                Fix(j,:) = Finew(j,l) ./ Fixsum(j);
                
                for k = 1:gl

                    if k == 1
                        
                        % Compute new cummulative grain size distributions
                        CummDist(j,k) = 0;

                    elseif k == 2

                        % Compute new cummulative grain size distributions
                        CummDist(j,k) = Fix(j,k-1);

                    else
                        
                        % Compute new cummulative grain size distributions
                        CummDist(j,k) = CummDist(j,k-1) + Fix(j,k-1);

                    end                   

                end

                CD1 = find(CummDist(j,:)<0.90,1,'last');         
                CD2 = find(CummDist(j,:)>0.90,1,'first');
                D90x(j) = (2 ^ (PsiScale(1,CD1) + (0.90 - CummDist(j,CD1)) * ((PsiScale(1,CD2) - PsiScale(1,CD1)) / (CummDist(j,CD2) - CummDist(j,CD1))))) / 1000;
                Lax(j) = D90x(j) * GrainScaleConstant;
                Ksx(j) = Nk * D90x(j);
                Fs1 = find(Grain==2,1,'last');
                Fsx(j) = CummDist(j,Fs1);
                % EXECUTE THE SURFACE MATERIAL GRAINSIZESTATISTICSWORKING1 FUNCTION
                Fipass(1,:) = Fix(j,:);
                [ Sigmag,Dg ] = GrainSizeStatistics( Gl,Psi,Fipass );   
                Dgx(j) = Dg;
                Dgmm(j) = Dgx(j) * 1000;
                DgSq1x(j) = Dgmm(j)^2;
                RMSDgmm = sqrt((Nn^-1.*sum(DgSq1x)));
                Sigmax(1,j) = Sigmag;

            end
            
            % Fill the grain size matrix with the new surface distributions
            % On the next time step the updated matrix will be directed to
            % teh Mobility analysis.  Keep the first two rows from the
            % original matrix to avoid changing the Mobility code.
            %GSD(:,3:end) = CummDist';
            
            % Compute active layer difference for use in surface material
            % expression in the next time step.  Use linear indexing.
            j =1:Nn;
            DLax = (Lax(j) - Laxp(j)) ./ Lax(j);
            
        end
         
        %% SAVE RESULTS OF BED EVOLUTION AT A FEW INSTANCES IN TIME FOR ANALYSIS
        if rem(i,timesave) == 0
            
            disp(i)
                
            time(index,1) = i * dt / 3600;
            Fisave(index,:,:) = Fix(:,:);
            Ssave(index,:) = So(1,:);
            BedSlope(index,:) = S(1,:);
            etasave(index,:) = n1(1,:);
            Dgsave(index,:) = Dgx(1,:);
            Sigmasave(index,:) = Sigmax(1,:);
            D90save(index,:) = D90x(1,:);
            WSEsave(index,:) = WSE(1,:);
            ysave(index,:) = y(1,:);
            Vsave(index,:) = V(1,:);
            Ushearsave(index,:) = Ushear(1,:);
            if SED == 2 || 3
                RMSSize(index,:) = RMSDgmm(1,:);
            end
            Qbsave(index,:) = Qbe(1,:);  
            Fssave(index,:) = Fsx(1,:);
            TotalTime = toc(TimeStart);
            if HYD == 2
                AdvCounterSave = AdvCounter;
            end 
            index = index + 1;
            
        end
        
        waitbar(i/tt)

    % Fill the n1 vector with the new bed elevations prior to next time iteration    
    n1 = etanew(1,:);

    end

    close(h)

    TotalTime = toc(TimeStart);

    %% WRITE DATA TO FILES FOR ANALYSIS
        
    save('stationing.mat','Station');
    save('time.mat','time');
    save('grainsizedistribution.mat','Fisave');
    save('slope.mat','Ssave');
    save('bedslope.mat','BedSlope');
    save('bedelevation.mat','etasave');
    save('grainsize.mat','Dgsave');
    save('standarddeviation.mat','Sigmasave');
    save('d90save.mat','D90save');
    save('watersurface.mat','WSEsave');
    save('depth.mat','ysave');
    save('Velocity.mat','Vsave');
    save('ushear.mat','Ushearsave');
    if SED == 2 || SED == 3
        save('rmsmean.mat','RMSSize');
    end
    save('transportrates.mat','Qbsave');
    save('sandfraction.mat','Fssave');
    save('elapsedtime.mat','TotalTime');
    if HYD == 2
        save('advcounter.mat','AdvCounterSave');
    end

end

