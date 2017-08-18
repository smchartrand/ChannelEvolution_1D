function [ TCMPercentileMatrix  ] = FAMobilityModelSimple( GSD )

    %% THIS CODE EXPLORES HOW friction ANGLE DISTRIBUTIONS VARY FOR
    % FOR DIFFERENT GRAIN SIZE DISTRIBUTIONS. THE CODE WAS MOTIVATED
    % KIRCHNER ET AL., 1990 AND BUFFINGTON ET AL., 1992.

        %% LOAD A BED SURFACE GRAIN SIZE DISTRIBUTION FOR FRICTION ANGLE 
        % analysis from the bed evolution model. The first column is the grian diameter 
        % in millimeter, the second column is grain diameter in Psi units, the third column 
        % is the cummulative distribution. Each spatial node is analyzed separately.
        % Grain size diameters in mm.
        GSDdiameter = GSD(:,1); %#ok<*NASGU>
        % Grain size diameters in psi units.
        GSDpsi = GSD(:,2);
        % All the cummulative distributions.
        gsd = GSD(:,3);
        % Specify the number of grains to use to build the virtual streambed
        VBed = 1000;
        VBedPlus = VBed+1;
        % Specify the number of random grains to place on the bed to measure
        % friction angles
        VRandom = 5000;
            
            % Pass the correct grain size cumulative distribution into
            % the calculations - line up grain size distribution with the
            % correct spatial node
            gsdpass = gsd;           
            
            %% IDENTIFY THE PsiDIFF TERM BASED ON THE SELECTED GRAIN SIZE DISTRIBUTION
            % First identify the address of the 'last' 0 value in the cumm
            % distribution of the selected grain size distribution and then
            % the 'first' 1 in the distribution.  This eliminates size
            % classes which do not account for mass within the distribution
            g0 = find(gsdpass==0,1,'last');
            g100 = find(gsdpass>0.9999,1,'first');
            % Changed this to operate on the array rather than row or
            % column for each spatial node.
            GSDuse = gsdpass(g0:g100);
            GSDpsiuse = GSDpsi(g0:g100);

            % Now use the results from the specified grain size distribution to
            % figure out some details needed for the calculations below.

            % Length of the grain size distribution vector
            GSDlength = length(GSDuse) - 1;
            % Specify the minimum grain size in psi units
            PsiMIN = min(GSDpsi (g0));
            % Specify the maximum grain size in psi units
            PsiMAX = max(GSDpsi(g100));
            % Find the difference between the max and min psi values
            PsiDIFF = PsiMAX - PsiMIN;

            %% COMPUTE SOME GRAIN SIZE DISTRIBUTION DETAILS
            % I have vectorized this operation for speed; the vector variable is k
            k = 1:GSDlength;

            % Compute the fractional percent of material in the bed
            Fi(k) = GSDuse(k+1) - GSDuse(k); 
            % Compute characteristic grain size of each size class in the bed
            Psi(k) = (GSDpsiuse(k) + GSDpsiuse(k+1)) * 0.5;

            % CALL THE GRAIN STATISTICS FUNCTION TO COMPUTE THE D50 AND D84
            % D50 and D84 are returned in meters
            [ D50,D84 ] = GrainSizeStatisticsMOB( GSDlength,Psi,Fi,GSDuse,GSDpsiuse );

            %% STEP ONE: POPULATE A LINE SEGMENT WITH VBed GRAINS.  
            % Grain size is determined by a random number, and grain coordinates 
            % are specified by their respective diameters. This is three step
            % process.

            % STEP ONE: generate a set of random numbers. Shuffle indicates that
            % the random number generator will generate a different sequence of
            % numbers each time rand is called.

                rng('shuffle');

                % Populate the R vector with VBedPlus random numbers. The R1 vector requires
                % VBedPlus grains so that a grain placement address exists just before the
                % last grain b/c I assume that grains will come to rest at the edge
                % between two adjacent grains. R1 is used to determine a random grain 
                % size.
                R1 = rand(VBedPlus,1);
                R1l = length(R1);

            %% STEP TWO: COMPUTE RANDOM GRAIN SIZES FROM THE RANDOM NUMBER SEQUENCE.
            % Use the random number array to specify random grain sizes to populate 
            % the computational array.  Grain size is in millimeter. 

                % R1GSet is the elevation of the R1GS grain top in millimeter.
                R1GSet = zeros(VBedPlus,1);

                % Compute R1Psi - R1Psi stands for random grain size in psi units.
                R1Psi = PsiMIN + (PsiDIFF .* R1);
                % Compute R1GS - this is equivalent to the inverse ln function.
                % Convert the random grain sizes from psi units to millimeter. 
                % R1GS stands for random grain size in millimeter.
                R1GS = 0.999775 * exp(0.693822.*R1Psi);

            %% STEP THREE: POPULATE THE LINE SEGMENT WITH THE RANDOM GRAIN SIZES.
            % Pre-populate the location and elevation vecotrs. centx is the spatial
            % coordinate of the grain center along the line segment. Center 
            % coordinate is in millimeter.

                centx = zeros(VBed,1);
                % edgex is the spatial coordiante of far grain edge along the line 
                % segment. Edge coordinate is in millimeter.
                edgex = zeros(VBed,1);
                % e0 is the elevation of each grain along the grain size vector, set 
                % arbitrarily to 1 so that the center of each grain lies in the same 
                % horizontal plane. Elevation is in millimeter.
                e0 = zeros(VBedPlus,1);

                MovingStaCumm = zeros(VBedPlus,VBedPlus);
                MovingStaCummRev = zeros(VBedPlus,VBedPlus);
                constant = zeros(VBedPlus,1);
                constantrev = zeros(VBedPlus,1);

                % Now run through and compute parameter values in order from 1 to
                % the end of the grain size vector
                for i = 1:VBedPlus

                   if i == 1
                       edgex(i) = 0 + R1GS(i);
                       centx(i) = 0 + (R1GS(i) / 2);
                       e0(i) = 1;
                       R1GSet(i) = e0(i) + (R1GS(i) / 2);

                   else
                       edgex(i) = edgex(i-1) + R1GS(i);
                       centx(i) = edgex(i-1) + (R1GS(i) / 2);
                       e0(i) = 1;
                       R1GSet(i) = e0(i) + (R1GS(i) / 2);

                   end

                end

            %% STEP FOUR: USING RESULTS OF STEP THREE, COMPUTE FORWARD MOVING WINDOW
            % matrices and backward moving window matrices of cummulative
            % stationing. this is needed to compute average elevation around any
            % grain based on Kirchner et al., 1990 because I will need to search
            % around a grain for a distance of D84 to compute the average grain
            % elevation.

                % I will read the edgex vector to a new term to keep naming
                % convention consistent through the next several steps. 
                % StaCumm is the cummulative stationing of the random grain size
                % vector.  StaCum is in millimeters.  This term will be used to
                % help compute the average grain top elevation for use in Kirchner
                % et al, 1990.
                StaCumm = edgex;

                % Now compute the moving window matrix of cummulative position
                for h = 1:VBedPlus

                   % This if statement creates a matrix of values which provides a
                   % moving window of cummulative stationing for the R1GS vector.
                   % Negative values in any vector of the matrix will be ignored in
                   % the average elevation calculations.
                   if h == 1

                       % Need this term to create the matrix
                       constant(h) = StaCumm(h);
                       % This operation subtracts the constant value from each 
                       % element of the MovingStaCumm vector.               
                       MovingStaCumm(h,:) = StaCumm - constant(h);

                   else

                       % Need this term to create the matrix
                       constant(h) = StaCumm(h);
                       % This operation subtracts the constant value from each 
                       % element of the MovingStaCumm vector.               
                       MovingStaCumm(h,:) = StaCumm - constant(h);

                   end  

                end

                % Now run through and do the above operations in reverse in order 
                % to compute elevations for Kirchner et al., 1990
                for i = VBedPlus:-1:1

                   if i == VBedPlus
                       edgexrev(i) = 0 + R1GS(i);

                   else
                       edgexrev(i) = edgexrev(i+1) + R1GS(i);

                   end

                end

                % I will read the edgexrev vector to a new term to keep naming
                % convenction consistent through the next several steps.
                StaCummRev = edgexrev;

                for h = VBedPlus:-1:1

                   % This if statement creates a matrix of values which provides a
                   % moving window of cummulative stationing for the R1GS vector.
                   % Negative values in any vector of the matrix will be ignored in
                   % the average elevation calculations.
                   if h == VBedPlus

                       % Need this term to create the matrix
                       constant(h) = StaCummRev(h);
                       % This operation subtracts the constant value from each 
                       % element of the MovingStaCumm vector.               
                       MovingStaCummRev(h,:) = StaCummRev - constant(h);

                   else

                       % Need this term to create the matrix
                       constant(h) = StaCummRev(h);
                       % This operation subtracts the constant value from each 
                       % element of the MovingStaCumm vector.               
                       MovingStaCummRev(h,:) = StaCummRev - constant(h);

                   end  

                end

            %% STEP FIVE: COMBINE THE FORWARD AND BACKWARD LOOKING MATRICES INTO 
            % one matrix for later calculations.  I will replace negative values in
            % the forward matrix with positve values from the backward matrix. This
            % is done with logical indexing; I replace the negative components of 
            % the forward looking matrix with positve values of the backward
            % looking matrix.

                MovingStaCumm(MovingStaCumm < 0) = MovingStaCummRev(MovingStaCummRev > 0);

            %% STEP SIX: NOW CREATE A VECTOR OF AVERAGE BED ELEVATION BASED ON A 
            % search neighborhood equivalent to the D84.

                % AvgElevR1GS is the average elevation of the R1 grain size vector
                % based on a search neighborhood of D84.  units are in millimeters.
                AvgElevR1GS = zeros(VBedPlus,1);

                for h = 1:VBedPlus

                    [~,c] = find(MovingStaCumm(h,:) < (D84*1000));

                    AvgElevR1GS(h) = sum(R1GSet(c)) / length(c);                 

                end

            %% PLACE GRAINS ALONG THE LINE SEGMENT AND COMPUTE friction ANGLES
            % Using the array of VBed random grains placed along the line segment,
            % randomly place grains of random size and compute the friction angle for
            % each occurrence. This is a three step process.

            % STEP ONE: generate two sets of random numbers. Shuffle indicates
            % that the random number generator will generate a different sequence 
            % of numbers each time rand is called.

                rng('shuffle');

                % R2 is used to determine a random grain size.
                R2 = rand(VRandom,1);
                % R3 is used to determine a random location along the grain size 
                % array populated with grains specified by R1. R3 is generated with
                % the random integer generator with replacement which means the 
                % same location can be generated more than once.  Note that there
                % is N-1 locations along the array specified above because it is 
                % not possible to compute a friction angle for the last grain in the
                % sequence because there is no grain down array from it.
                R3 = randi([1,VBed],[VRandom,1]);
                % R4 is used to determine a new grain size to place along the grain
                % size array for determination of a friction angle.  R4 is rounded to
                % the nearest integer value and mapped to spatial location 1, 2, 3,
                % etc. along the grain size array. 
                R4 = rand(VRandom,1);

            % STEP TWO: locate the R2 random grain sizes along the R1 array
            % according to the random location specified by R3.  Then compute the
            % R3 grain elevation and finally the friction angle.  Write the friction
            % angle to a parameter and save.

                edgexlessone = zeros(VRandom,1);

                % R2Psi stands for random grain size in psi units.
                R2Psi = PsiMIN + (PsiDIFF .* R2);
                % Convert the random grain sizes from psi units to millimeter. R2GS stands 
                % for the randomly placed grain, size expressed in millimeter.
                R2GS = 0.999775 * exp(0.693822.*R2Psi);

            % STEP THREE: map the grain sizes from the R1 vector to the R3 location
            % vector for use in looping below and to determine where the placed 
            % grain can physically fit. To determine if the placed grain can
            % physically fit I need to store grain sizes beyond the two bounding
            % grains of the placement address because their geometry maybe
            % incompatible with the size of the placed grain.  In that case I need
            % to look down array to the next grain size, and perhaps beyond to
            % determine when the grain may physically fit.

                % edgexi term is the edgex location from the R1 vector mapped to
                % the R3 location. units are in millimeter.
                edgexi = edgex(R3);
                edgexii = edgex(R3+1);
                % centx2 term is the centx location from the R1 vector mapped to  
                % the R3 location, and so on. units are in millimeter.
                centxi = centx(R3);
                centxii = centx(R3+1);
                % R3i term is the R1GS grain size from the R1 vector mapped to the
                % R3 location.  This mapping is done by linear indexing with R3.
                % units are in millimeter.
                R3GSi = R1GS(R3);
                R3GSiet = (R3GSi ./ 2) + e0(1);
                % R3ii term is the R1GS+1 grain size from the R1 vector mapped to
                % the R3 location. This mapping is done by linear indexing with R3.
                % units are in millimeter.
                R3GSii = R1GS(R3+1);
                % R3AvgElevR1GS is the R1GS average elevation mapped to the R3
                % vector for use in comparing against the R2GS elevation.
                R3AvgElevR1GS = AvgElevR1GS(R3);
                % Run through a loop to map the grain sizes upstream of i in order
                % to compute the exposure term properly. the notation of putting i
                % before the parameter name is to indiacte one value up vector from
                % the i position.
                iR3GSi = zeros(VRandom,1);

                for i = 1:VRandom;

                    if R3(i) > 1

                        iR3GSi(i) = R1GS(R3(i)-1);   

                    end

                end

            % STEP FOUR: now roll through and compute friction angles and several 
            % other terms which are needed for evaluating the critical shear stress
            % function of Wiberg and Smith, 1987 and Kirchner et al, 1990.

                dl = zeros(VRandom,1);
                % dl is the length between the
                hyp = zeros(VRandom,1);
                % theta is the computed friction angle in degrees
                theta = zeros(VRandom,1);
                % R2GSe is the elevation of the R2 grain center placed at location
                % R3. R2GSet is the elevation of the R2 grain top.
                R2GSec = zeros(VRandom,1);
                R2GSet = zeros(VRandom,1);
                % only use the fist value in the e1 vector because e1 is 
                % constant. if this changes must deal with e1.
                E1(1) = e0(1);
                % exposure is the grain exposure for use in Kirchner1990.
                exposure = zeros(VRandom,1);
                % exposureelevdiff is the relative elevation difference between the
                % random R2 grain top and the upgrain R3i grain top
                exposureelevdiff = zeros(VRandom,1);
                % RelGS is the ratio of the R2 grain size to the D50 of the
                % mixture
                RelGS = zeros(VRandom,1);

                for i = 1:VRandom

                    % This if statement tests whether the upstream grain is 1.5x
                    % the next downstream grain, and whether the random
                    % grain to be placed is larger than the diameters of the u/s
                    % and d/s grains along the R vector.
                    if (R3GSi(i) / R3GSii(i)) >= 100000

                        if R2GS(i) < ((R3GSi(i) + R3GSii(i)) / 2)

                            % dl computation assumes the center of the overlying
                            % grain falls precisely at the edge location between 
                            % the R3i and R3ii grains.
                            dl(i) = (R3GSii(i,1) / 2);

                            % the hypotenuse is computed as the simple sum of the 
                            % radius's for the overlying and the R3ii grains.
                            hyp(i) = (R2GS(i) / 2) + (R3GSii(i,1) / 2);

                            % the angle theta is the friction angle.
                            theta(i) = asind(dl(i) / hyp(i));

                        else

                            theta(i) = 1;

                        end

                    elseif (R3GSii(i) / R3GSi(i)) >= 100000

                        if R2GS(i) < ((R3GSi(i) + R3GSii(i)) / 2)

                            % dl computation assumes the center of the overlying
                            % grain falls precisely at the edge location between 
                            % the R3i and R3ii grains.
                            dl(i) = (R3GSii(i,1) / 2);

                            % the hypotenuse is computed as the simple sum of the 
                            % radius's for the overlying and the R3ii grains.
                            hyp(i) = (R2GS(i) / 2) + (R3GSii(i,1) / 2);

                            % the angle theta is the friction angle.
                            theta(i) = asind(dl(i) / hyp(i));

                        else

                            theta(i) = 89;

                        end

                    else

                        % dl computation assumes the center of the overlying
                        % grain falls precisely at the edge location between 
                        % the R3i and R3ii grains. units are in millimeter.
                        dl(i) = (R3GSii(i,1) / 2);

                        % the hypotenuse is computed as the simple sum of the 
                        % radius's for the overlying and the R3ii grains.
                        % units are in millimeter.
                        hyp(i) = (R2GS(i) / 2) + (R3GSii(i,1) / 2);

                        % the angle theta is the friction angle.
                        theta(i) = asind(dl(i) / hyp(i));

                        % R2GSec is the elevation of the R2 grain center placed at 
                        % location R3 or R3new. R2GSet is the elevation of the R2 grain 
                        % top. units are in millimeter.
                        R2GSec(i) = ((cosd(theta(i)) * (R2GS(i) / 2)) + (sind(theta(i)) * (R3GSii(i) / 2))) + E1;
                        R2GSet(i) = R2GSec(i) + (R2GS(i) / 2);   

                        % Compute the relative grain size for R2GS
                        RelGS(i) = R2GS(i) ./ (D50 .* 1000);

                        % Compute the R2GS grain exposure for use in the Kirchner1990.
                        % units are in millimeter.  Exposure equals the top elevation
                        % of the random placed grain minus the top elevation of bed
                        % grain i.  This is a simplification of reality.  I should
                        % additionally search up grain from grain i to assess if grains
                        % extend above the random placed grain.  I will do that in
                        % the future.

                        exposureelevdiff(i) = R2GSet(i) - R3GSiet(i);

                        if exposureelevdiff(i) > 0

                            exposure(i) = R2GSet(i) - R3GSiet(i);

                        else

                            exposure(i) = 0;

                        end

                    end

                end

            %% CALL THE KIRCHNER FUNCTION TO COMPUTE CRITICAL DIMENSIONLESS 
            % shear stress distributions from the friction angle distributions

            % remember that some terms are in units of millimeters. convert to
            % meters inside the Kirchner funtion.

            % zo is the length scale of logarithmic velocity profile.  
            % units are in millimeters.  I will specify that zo is a function of the
            % local average bed elevation rather than Dietrich and Whiting,
            % 1989 because local conditions can vary significantly from the
            % static case
            zo = R3AvgElevR1GS;
            % set z equal to R2GSet. units are in millimeter.
            z = R2GSet;
            [taucrit,taucritless,R2GSeval,ycrit,Scrit] = Kirchner1990v3(R2GS,theta,exposure,z,zo,R3GSiet,VRandom);

            %% CALL THE MOBILITY FUNCTION TO COMPUTE PDFS OF MOBILITY
            % based on 0.5 PSI class intervals.
            [TCMPercentileA] = MobilityPDFMOB(GSDlength,taucrit,taucritless,PsiMIN,PsiMAX,R2GSeval,VRandom);
           
            % Build a matrix to store the critical dimensionless stress
            % percentile results by spatial node.  The results are stored
            % as row = spatial node and column = critical dimensionless
            % stress to pass to the sediment transport computations
            TCMPercentileMatrix = TCMPercentileA;

end

