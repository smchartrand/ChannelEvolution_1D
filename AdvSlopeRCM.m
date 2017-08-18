function [ yout,counter2 ] = AdvSlopeRCM( yout,j,counter2,yr,dx,HLimd,Smp )
    
    %% REDUCED COMPLEXITY APPROXIMATION OF ADVERSE SLOPE WATER DEPTH. 
    % This code employs a reduced complexity approximation of water depth within
    % pool segments that contain horizontal of adverse bed slopes.  The 
    % approximation is based on projecting an average bed slope over the reach
    % in order to compute a residual pool depth.  From there the approximation
    % computes the depth based on a slope factor which spans the range from 
    % zero to the reach average slope. The approximate water depth increases as 
    % the normal flow depth for the average bed slope increases towards a factor
    % two of the residual pool depth. The factor two represents the reduced 
    % complexity rule. The rule is based on the physical observation that
    % water surface slopes within pools at low flow are quite flat,
    % increasing as the discharge increases.
    
    % Initialize a second counter used in the adverse slope reduced
    % complexity model water depth function. The counter aids in 
    % computing the length upstream from the first node of 
    % adverse/horizontal slope. This ultimately leads to computation 
    % of the depth greater than the residual depth, defined by the 
    % projection of the average bed slope at the node of interest.
    if counter2 == 0
        
        counter2 = 1;
        
    else
        
        counter2 = counter2 + 1;
        
    end
    %% THE CALCULATION IS SIMPLE AND DEPENDS ON THE CHOSEN MINIMUM LIMITING POOL 
    % slope and how far upstream any particular node of negative or
    % horizontal slope is located from the downstream critical control.
        
    yout(1,j) = ((counter2 * dx) * tan(Smp))+ HLimd(1,j) + yr(1,j);
               
end

