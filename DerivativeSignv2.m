function [ order,nc,ncindex ] = DerivativeSignv2( j,n1,Nn,dx,order,nc,ncindex )

    %% UNTITLED Summary of this function goes here
    % Detailed explanation goes here
    
    %% STEP 1: FIND HIGH POINT LOCATIONS WITHIN THE BED ELEVATION VECTOR
    if j > 1 && j < Nn    
        
        % Flag the sign of the bed elevation derivative
        if ((n1(1,j-1) - n1(1,j)) / dx) < 0 && ((n1(1,j) - n1(1,j+1)) / dx) > 0
            
            order(1,j) = 1;
            nc(1,j) = n1(1,j);
            
        end
            
    end
    
    if j == Nn

        % Now flag all locations along the vector which are high points
        ncindex = find(nc ~= 0);
        L = length(ncindex);

        if L > 1

            for i = 1:L

                if i > 1

                    if nc(1,ncindex(i-1)) < nc(1,ncindex(i))

                        nc(1,ncindex(i-1)) = nc(1,ncindex(i));

                    end

                end


            end

        end

    end
    
end

