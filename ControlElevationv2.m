function [ yr ] = ControlElevationv2( n1,Nn,nc,ncindex,yr,S )

    %% UNTITLED Summary of this function goes here
    % Detailed explanation goes here
    
    l = length(ncindex);
    
    %% STEP 1: FIND HIGH POINT LOCATIONS WITHIN THE BED ELEVATION VECTOR
    for i = 1:l
        
        for j = 1:Nn
        
            if i == 1

                if j < ncindex(i) && n1(1,j) < nc(1,ncindex(i))

                    nc(1,j) = nc(1,ncindex(i));
                    
                end

            elseif i > 1

                if j > ncindex(i-1) && j < ncindex(i) && n1(1,j) < nc(1,ncindex(i))

                    nc(1,j) = nc(1,ncindex(i));


                end

            end
            
            % Compute the residual water depth (yr) - meters
            % The residual water depth is used in the RCM calculation of
            % pool water depth along segments of horizontal or negative bed
            % slope.
            if S(1,j) < 0 && nc(1,j) ~= 0
            
                yr(1,j) = nc(1,j) - n1(1,j);
                
            end
        
        end
        
    end
    
end

