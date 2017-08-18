function [ Qbe ] = AshidaMichiuev2( Nn,Gl,Fix,Dgx,Ushear,GSm  )

    %% ASHIDA AND MICHIUE FUNCTION COMPUTES BEDLOAD TRANSPORT
    % This module computes fractional bedload transport using the Ashida-
    % Michiue equation. The script has been vectorized to optimize speed as 
    % it is computationally intensive for large ranges of grain size.  See 
    % Ashida-Michiue (1972) Transactions, Japan Society for Civil Engineering
    % for details.  Formulation of the Ashida-Michiue function is based on the Parker 
    % e-book Chapter 7, which differs from the published formulation in 
    % that the Parker version uses the geometric mean of the surface grain size,
    % as opposed to the arithmatic mean.

    %% PRE-ALLOCATE MEMORY FOR PARAMETERS    
    PBKtranspose(Nn,Gl)=0;
    Qbe(1,Nn)=0;
    tstarrkc(1,Gl)=0;
    tstarrk(1,Gl)=0;
    Dratio(1,Gl)=0;
    Dratioo(1,Gl)=0;
    Dratiooo(1,Gl)=0;

    %% DEFINE A FEW CONSTANTS RATHER THAN PASSING THEM TO THE FUNTION
    % Relative density - dimensionless
    Rd = (2650/1000) - 1;   

    % Gravitational acceleration - meters per square second
    g = 9.81;      

    %% COMPUTE THE SEDIMENT TRANSPORT RATES BY PARTICLE SIZE
    % j is the spatial node
    for j = 1:Nn
    
        % Pass values for each spatial node rather than indexing the 
        % calculations to each spatial node
        Ushearpass(1) = Ushear(1,j); 
        Fipass(1,:) = Fix(j,:);
        Dgpass(1) = Dgx(j);
        
        % Set the vector index parameter - the number of grain classes
        k = 1:Gl;
        
        % Compute the grain size specific Shields number
        Fipass1 = Fipass > 0;
        tstarrk(Fipass1) = (Ushearpass ^ 2) ./ (Rd .* g .* GSm(Fipass1));
        
        % The effective critical Shields stress for surface geometric size
        % according to Ashida and Michiue, 1972
        tstarrg = 0.05;

        % Compute the grain size ratio for transport computation
        Dratio(Fipass1) = GSm(Fipass1) ./ Dgpass;
        % Create an index of phi values < 1.35 for computing Wstark
        Dratioo(Fipass1) = Dratio(Fipass1) <= 0.40;
        ind1 = logical(Dratioo);
        Dratiooo(Fipass1) = Dratio(Fipass1) > 0.40;
        ind2 = logical(Dratiooo);
        
        % Compute dimensionless transport rate for indexed ind1 = true
        tstarrkc(ind1)= (0.843.*(Dratio(ind1).^ -1)) .* tstarrg;
        % Compute dimensionless transport rate for indexed ind2 = true
        tstarrkc(ind2) = (log10(19) ./ (log10(19 .* Dratio(ind2)))).^2 .* tstarrg;
        % Set dimensionless transport rate for all other size classes to
        % zero
        tstarrkc(~Fipass1) = 0;

        % Compute the fractional bedload transport by grain class -
        % dimensional - square meters per second
        qbstark = 17 .* (tstarrk - tstarrkc) .* (tstarrk .^ 0.5 - tstarrkc .^ 0.5);
        %qbstark(1,Fipass1) = 17 .* (tstarrk(Fipass1) - tstarrkc(Fipass1)) .* (tstarrk(Fipass1) .^ 0.5 - tstarrkc(Fipass1) .^ 0.5);
        
        % Compute the dimensional sediment transport by grain class (square
        % meters pers second)
        qbk = qbstark .* ((Rd .* g .* GSm(k)) .^ 0.5 .* GSm(k) .* Fipass);
        %qbk = qbstark(1,Fipass1) .* ((Rd .* g .* Dipass(Fipass1)) .^ 0.5 .* Dipass(Fipass1) .* Fipass(Fipass1));

        % Compute the total sediment transport for a specific node - 
        % dimensional - square meters per second
        qb = sum(qbk);

        % Compute the percent fractional bedload transport by grain class
        pbk = qbk ./ qb;

        % Store the total bedload transport by spatial node for passing
        Qbe(1,j) = qb;

        % Store the fractional bedload transport by grain class for passing
        PBKtranspose(j,:) = pbk; 
    
    end
    
end

