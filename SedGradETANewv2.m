function [ sedgrad,etanew] = SedGradETANewv2(alpha,Qbe,Qsf,dt,dx,Lamda,n1,No,Nn,HYD,DWSE,y)
    
    %% SEDIMENT TRANSPORT GRADIENT AND CHANNEL BED PROFILE EVOLUTION
    % This function computes the spatial sediment transport graident based
    % on central differencing and uses these results to compute channel bed
    % profile evolution along a 1-D channel.  The function has been
    % optimized for speed by indexing the spatial loop.  Indexing the
    % spatial loop is tricky and should be reviewed carefully for errors.
    
    % Set the vector index parameters - the number of spatial nodes
    j = 1:Nn;
    k = 1:Nn-1;
    l = 1:Nn-2;
    % Create a set of index terms to specify spatial location.  There are
    % many index terms so it is easy to get lost and make a mistake.  Step
    % through this part carefully.  Basically lining up vectors of data
    % with certain ranges of spatial nodes then truncating the vectors to
    % be all the same length in order to perform operations.
    % Just the first spatial node
    spatialindex1 = j == 1;
    % All nodes except the last spatial node
    spatialindex2 = j < Nn;
    % All nodes except the first node
    spatialindex3 = j > 1 & j <= Nn;
    % All nodes except the first and last node
    spatialindex4 = j > 1 & j < Nn;
    % Just the last spatial node
    spatialindex5 = j == Nn;
  
    % Now extract the sediment transport rates that line up with all nodes
    % except the last spatial node
    Qbeindex2 = Qbe(spatialindex2);
    % Now do the same for the sediment transport rates.  This gives me a
    % vector of transport rates for cells j = 1 through Nn - 1.
    Qbeindex2(Qbeindex2 == 0) = [];
    % Now extract the sediment transport rates that line up with all nodes
    % except the first spatial node
    Qbeindex3 = Qbe(spatialindex3);
    % Now do the same for the sediment transport rates.  This gives me a
    % vector of transport rates for cells j = 2 through Nn.
    Qbeindex3(Qbeindex3 == 0) = [];
    
    % Now compute the sediment transport gradient for all nodes.  The
    % gradient is defined as the difference from downstream to upstream.
    % For example, the gradient operating on the first two nodes is
    % computed as node 2 (downstream) less node 1 (upstream).
    Qbegrad1(k) = Qbeindex3(k) - Qbeindex2(k);
    % Now delete the first and last cells of the sedgrad1 array in order
    % to develop arrays that can be operated on within a vectorized
    % statement.  The first array will contain sedgrad values for nodes 2
    % through Nn - 1.  The second array will contain sedgrad values for
    % nodes 3 through Nn.  Node 2 equals node 2 less node 1 and Node 3
    % equals node 3 less node 2.  First make 2 copies of the sedgrad1
    % array.
    Qbegrad11 = Qbegrad1;
    Qbegrad111 = Qbegrad1;
    % Now clip the last and first elements of the sedgrad1 array for use in
    % the operations below.  Now I have all the data needed to compute the
    % sediment gradients.
    Qbegrad11(end) = [];
    Qbegrad111(1) = [];
    
    % Now move through an compute the sediemnt gradients.  Sediment transport 
    % gradient for the first spatial node based on the boundary condition.  This 
    % is the left hand term of the sediment gradient equation.
    sedgrad1(spatialindex1) = (alpha .* (Qbe(spatialindex1) - Qsf));
    % Sediment transport gradient for off set nodes + 1.  This is the right
    % hand term of the sediment gradient equation
    sedgrad11(k) = (1-alpha) .* (Qbegrad1);
    % Sediment transport gradient at the first spatial node
    sedgradN1 = sedgrad1 + sedgrad11(1);
    
    % Sediment trasnport gradient for nodes 2 through Nn - 1.  First deal
    % with the left hand term and second deal with the right hand term.
    % Left hand term
    sedgrad111(l) = (alpha) .* (Qbegrad11(l));
    % Right hand term
    sedgrad1111(l) = (1-alpha) .* (Qbegrad111(l));
    % Sediment transport gradient at spatial nodes 2 through Nn - 1.
    sedgradN2 = sedgrad111 + sedgrad1111;
    
    % Sediment transport gradients for node Nn
    sedgradN3 = (alpha .* (Qbegrad1(end)));
    
    % Now build one array with all the sedgradN values.  Array built in
    % order from N1 to N3.  First set the size of the sedgrad array to fill
    % it properly and manage its memory allocation
    sedgrad(1,Nn)=0;
    sedgrad(spatialindex1) = sedgradN1;
    sedgrad(spatialindex4) = sedgradN2;
    sedgrad(spatialindex5) = sedgradN3;
    % Now use the sedgrad array to build the etanew array.  First set the size of 
    % the etanew array to fill it properly and manage its memory allocation
    etanew(1,Nn)=0;
    % First compute all the new bed elevations.  Next I will replace the
    % last bed elevation value depending upon the downstream boundary
    % condition.
    etanew(j) = n1(1,j) - ((1 / (1 - Lamda)) * (dt / dx) * sedgrad(1,j));   
    % Now replace the last value of the etanew array depending on the
    % downstream boundary condition
    if HYD == 1 
        etanew(spatialindex5) = No;
    % If backwater, bed elevation at Nn is a function of sediment transport gradient
    elseif HYD == 2 || 3
        etanew(spatialindex5) = No;%DWSE - y(1,Nn);%n1(1,spatialindex5) - ((1 / (1 - Lamda)) * (dt / dx) * sedgrad(1,spatialindex5));
    end                  

end

