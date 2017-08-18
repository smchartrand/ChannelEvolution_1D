function [ Qbe,PBKtranspose ] = WilcockCrowe2( Nn,Gl,Fix,Fsx,Psi,Dgx,Tbede,Ushear,Density,B )

    %% WILCOCK AND CROWE FUNCTION COMPUTES BEDLOAD TRANSPORT
    % This module computes fractional bedload transport using the Wilcock-
    % Crowe equation. The script has been vectorized to optimize speed as 
    % it is computationally intensive for large ranges of grain size.  See 
    % Wilcock-Crowe (2003) Journal of Hydraulic Engineering for details.  
    % Formulation of the Wilcock-Crowe function is based on the Parker 
    % e-book Chapter 17, which differs from the published formulation in 
    % that the Parker version uses dimensionless values of the pertinent 
    % parameters, the published version uses dimensional values.

    %% PRE-ALLOCATE MEMORY FOR PARAMETERS
    Wstark(Gl)=0;
    Wstarksave(1,Nn)=0;
    PBKtranspose(Nn,Gl)=0;
    Qbe(1,Nn)=0;
    QbeWC(1,Nn)=0;
    Fspass(1)=0;
    Dgpass(1)=0;

    %% DEFINE A FEW CONSTANTS RATHER THAN PASSING THEM TO THE FUNTION
    % Relative density - dimensionless
    Rd = (Density/1000) - 1;   

    % Gravitational acceleration - meters per square second
    g = 9.81;      

    %% COMPUTE THE SEDIMENT TRANSPORT RATES BY PARTICLE SIZE
    % j is the spatial node
    for j = 1:Nn
    
        % Pass values for each spatial node rather than indexing the 
        % calculations to each spatial node
        Tbedepass(1) = Tbede(1,j);
        Ushearpass(1) = Ushear(1,j); 
        Fipass(1,:) = Fix(j,:);
        Fspass(1) = Fsx(j);
        Dgpass(1) = Dgx(j);

        % Set the vector index parameter - the number of grain classes
        k = 1:Gl;
        
        % Compute the reference shear stress for the surface Dg - dimensionless
        tstarrg = 0.021 + (0.015 * exp(-20 * Fspass));       

        % Compute the hiding function power by grain class - dimensionless
        b = 0.67./(1 + exp(1.5-(((2.^ Psi(1,k))./1000)./ Dgpass)));

        % Compute the reference stress ratio by grain class - dimensionless
        trefdimencrit = (((2.^ Psi(1,k))./1000)./ Dgpass).^ -b(1,k);
        
        % Compute the reference stress for grain class i - dimensionless
%         tstari = tstarrg .* trefdimencrit(1,k);

        % Compute the transport stress ratio for transport computation
        phi = (Tbedepass./tstarrg) .* trefdimencrit(1,k);
%         phi = (Tbedepass./tstari);
        %phi = (Tbedepass./tstarrg).*(trefdimencrit(1,k));
        % Create an index of phi values < 1.35 for computing Wstark
        phii = phi < 1.35;
        
        % Compute dimensionless transport rate for indexed phi = true
        Wstark(phii) = 0.002.*(phi(phii).^ 7.5);
        % Compute dimensionless transport rate for indexed phi = false
        Wstark(~phii) = 14.*((1 - (0.894./(phi(~phii).^ 0.5))).^ 4.5);
        
        Wstarksave(1,j) = sum(Wstark);

        % Compute the fractional bedload transport by grain class -
        % dimensional - square meters per second
        qbk = (Fipass(1,k).*Wstark(1,k).*(Ushearpass ^ 3))./(Rd .* g);

        % Compute the total sediment transport for a specific node - 
        % dimensional - square meters per second
        qb = sum(qbk);

        % Compute the percent fractional bedload transport by grain class
        pbk = qbk ./ qb;

        % Store the total bedload transport by spatial node for passing
        Qbe(1,j) = qb;
        
        % Convert transport rate to kilograms per minute
        QbeWC(1,j) = Qbe(1,j) .* B(1,j) .* 60 .* Density;

        % Store the fractional bedload transport by grain class for passing
        PBKtranspose(j,:) = pbk; 
    
    end
    
end

