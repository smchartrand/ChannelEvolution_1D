function [ D50,D84 ] = GrainSizeStatisticsMOB( GSDlength,Psi,Fi,GSDuse,GSDpsiuse )

%% GRAIN SIZE STATISTICS FUNCTION
% GrainSizeStatistic computes various statistics for a grain size 
% ditribution. Results from these computations are used in a mixed grain 
% size sediment transport model.  The script has been optimized for speed
% vectorizing all operations.

    % Set the vector index parameter - the number of grain classes    
    k = 1:GSDlength;

    % Compute average grain size for class expressed - psi units
    PsiAvg1(k) = Psi(k) .* Fi(k);

    % Compute the average grain size expressed - psi units
    PsiAvg = sum(PsiAvg1);
    
    % Compute the D50 grain size of the distribution - meters
    CD1 = find(GSDuse<0.50,1,'last');
    CD2 = find(GSDuse>0.50,1,'first');
    D50 = (2 ^ (GSDpsiuse(CD1,1) + (0.50 - GSDuse(CD1,1)) * ((GSDpsiuse(CD2,1) - GSDpsiuse(CD1,1)) / (GSDuse(CD2,1) - GSDuse(CD1,1))))) / 1000;
    
    % Compute the D84 grain size of the distribution - meters
    CD1 = find(GSDuse<0.84,1,'last');
    CD2 = find(GSDuse>0.84,1,'first');
    D84 = (2 ^ (GSDpsiuse(CD1,1) + (0.84 - GSDuse(CD1,1)) * ((GSDpsiuse(CD2,1) - GSDpsiuse(CD1,1)) / (GSDuse(CD2,1) - GSDuse(CD1,1))))) / 1000;
 
end

