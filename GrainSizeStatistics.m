function [ Sigmag,Dg ] = GrainSizeStatistics( Gl,Psi,Fipass )

%% GRAIN SIZE STATISTICS FUNCTION
% GrainSizeStatistic computes various statistics for a grain size 
% ditribution. Results from these computations are used in a mixed grain 
% size sediment transport model.  The script has been optimized for speed
% vectorizing all operations.

    % Set the vector index parameter - the number of grain classes    
    k = 1:Gl;

    % Compute average grain size for class expressed - psi units
    PsiAvg1(k) = Psi(k).*Fipass(k);

    % Compute the average grain size expressed - psi units
    PsiAvg = sum(PsiAvg1);

    % Compute the standard deviation for each class - Psi units
    Sigma1 = (Psi(k) - PsiAvg).^2 .* Fipass(k);  
    
    % Compute the standard deviation of the grain size distribution - Psi units
    Sigma = sum(Sigma1).^(1/2);

    % Compute the geometric standard deviation of the grain size
    % distribution - meters
    Sigmag = (2.^Sigma)./1000;
    
    % Compute the geometric mean grain size of the distribution - meters
    Dg = (2.^PsiAvg)./1000;
 
end

