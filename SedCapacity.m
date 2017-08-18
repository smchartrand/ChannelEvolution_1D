function [ Qbe ] = SedCapacity( Nn,Dgx,Tbede,Tcrit,B,Density )
    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION

    % Relative density - dimensionless
    Rd = (2650/1000) - 1;

    % Gravitational acceleration - meters per square second
    g = 9.81; 
    
    %% PRE-ALLOCATE LOOPING PARAMETERS AND OTHER PARAMETERS
    
    Qbstar(1,Nn)=0;
    Qbe(1,Nn)=0;

    %% COMPUTE THE RELEVANT HYDRAULICS

    % Compute node to node bed slope using forward, backward and central differences

    for j = 1:Nn
        
        Tbedepass(1) = Tbede(1,j);
    
        % Compute sediment transport capacity - dimensionless
        if Tbedepass - Tcrit <= 0

            Qbstar(1,j) = 0;

        else

            Qbstar(1,j) = 3.97 * (Tbedepass - Tcrit)^1.5;

        end

        % Compute sediment transport capacity - dimensional 
        Qbe(1,j) = Qbstar(1,j) * (((Rd)* g * Dgx(1,j))^0.5 * Dgx(1,j));
        
        % Convert transport rate to kilograms per minute
        QbeES(1,j) = Qbe(1,j) .* B(1,j) .* 60 .* Density;
    
    end

end

