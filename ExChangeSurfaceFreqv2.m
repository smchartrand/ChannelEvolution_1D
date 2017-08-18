function [ Finew,Laxp ] = ExChangeSurfaceFreqv2( sedgrad,PBKtranspose,pLi,n1,fi,Fix,Qbe,Qsf,dx,dt,Lax,DLax,i,Nn,alpha,alphased,Gl,etanew,Lamda )

    %% THIS FUNCTION COMPUTES THE EXCHANGE FREQUENCY AND NEW GRAIN SIZE
    % distribution in the bed based on the exchange frequency parameter, the thickness 
    % of the active layer and the fractional sediment transport gradient.

    fli(Nn,Gl)=0;
    flii(1,Gl)=0;
    fracsedgrad(Nn,Gl)=0;
    Finew(Nn,Gl)=0;

    for j = 1:Nn

        %% COMPUTE THE EXCHANGE FRACTIONS
        % Assign passing varaibles in order to have the indexing below
        % work.  This avoids an additional for loop.
        etanewpass = etanew(j);
        n1pass = n1(j);
        fipass = fi(j,:);
        Fixpass = Fix(j,:);
        PBKtransposepass = PBKtranspose(j,:);
        % Set the vector index parameter - the number of grain classes
        k = 1:Gl;
        % Compute the difference between the new bed elevations and the
        % old bed elevations one element per loop.  This is the indexing
        % parameter.
        etadiff = etanewpass - n1pass;

        % Create an index term for the exchange frequency.
        etadiffindex = etadiff < 0;
        % Pass the index term result to an array with Gl number of
        % elements.
        etadiffindexx = repmat(etadiffindex,1,Gl);

        % Assign the bed material exchange frequency parameter based on
        % whether the bed is degrading or aggrading.
        if etadiffindex == 1
            
            % Specify the bed material for index = true
            flii(etadiffindexx) = fipass(k);
            
        else

            % Specify the bed material for index = false
            flii(~etadiffindexx) = alphased .* Fixpass(k) + ((1 - alphased) .* PBKtransposepass(k));

        end
            
        % Now fill the fli matrix with the updated array results
        fli(j,:) = flii(:);

        %% COMPUTE THE NEW SURFACE FREQUENCIEs        
        % Compute the fractional sediment transport gradient using the
        % index     
        
        if j == 1

            fracsedgrad(j,k) = ((alpha./dx) .* ((Qbe(1,j) .* PBKtranspose(j,k)) - (Qsf .* pLi(k)))) + (((1-alpha)./dx) * ((Qbe(1,j+1) .* PBKtranspose(j+1,k)) - (Qbe(1,j) .* PBKtranspose(j,k))));

            if i == 1

                Finew(j,k) = Fix(j,k) + ((dt ./ Lax(1,j)) .* ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - (Fix(j,k) - fli(j,k));

            else

                Finew(j,k) = Fix(j,k) + ((dt ./ Lax(1,j)) .* ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - ((Fix(j,k) - fli(j,k)) .* DLax(1,j));

            end

            spatialindex = Finew(j,k) < 0;
            
            Finew(j,spatialindex) = 0;

        elseif j > 1 && j <= Nn - 1

            fracsedgrad(j,k) = ((alpha./dx) .* ((Qbe(1,j) .* PBKtranspose(j,k)) - (Qbe(1,j-1) .* PBKtranspose(j-1,k)))) + (((1-alpha)./dx) .* ((Qbe(1,j+1) .* PBKtranspose(j+1,k)) - (Qbe(1,j) .* PBKtranspose(j,k))));

            if i == 1

                Finew(j,k) = Fix(j,k) + ((dt ./ Lax(1,j)) .* ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - (Fix(j,k) - fli(j,k));

            else

                Finew(j,k) = Fix(j,k) + ((dt / Lax(1,j)) * ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - ((Fix(j,k) - fli(j,k)) .* DLax(1,j));

            end
            
            spatialindex = Finew(j,k) < 0;
            
            Finew(j,spatialindex) = 0;

        else

            fracsedgrad(j,k) = (alpha./dx) .* ((Qbe(1,j) .* PBKtranspose(j,k)) - (Qbe(1,j-1) .* PBKtranspose(j-1,k)));

            if i == 1

                Finew(j,k) = Fix(j,k) + ((dt ./ Lax(1,j)) .* ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - (Fix(j,k) - fli(j,k));

            else

                Finew(j,k) = Fix(j,k) + ((dt ./ Lax(1,j)) .* ((1 ./ (1 - Lamda)) .* (-fracsedgrad(j,k) + (fli(j,k) .* sedgrad(1,j))))) - ((Fix(j,k) - fli(j,k)) .* DLax(1,j));

            end

            spatialindex = Finew(j,k) < 0;
            
            Finew(j,spatialindex) = 0;

        end

    end

    % Save current active layer vector to compute active layer term
    % of surface material expression in the next time step
    Laxp = Lax;

end

