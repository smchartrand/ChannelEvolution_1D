function [ Tbede,Ushear,y,WSE,V ] = Normal( n1,dx,Ksx,Dgx,Nn,qw )
    %% DEFINE SOME CONSTANTS RATHER THAN READING THEM IN TO THE FUNCTION

    % Relative density - dimensionless
    Rd = (2650/1000) - 1;             

    % Water density - kg per cubic meter
    Rho = 1000;

    % Gravitational acceleration - meters per square second
    g = 9.81; 
    
    % Manning-Strickler alpha parameter - dimensionaless
    a = 8.1;
    
    %% PRE-ALLOCATE LOOPING PARAMETERS AND OTHER PARAMETERS
    S(1,Nn)=0;
    y(1,Nn)=0;  
    WSE(1,Nn)=0;
    V(1,Nn)=0;
    Tbed(1,Nn)=0;
    Tbede(1,Nn)=0;
    Ushear(1,Nn)=0;

    %% COMPUTE THE RELEVANT HYDRAULICS

    % Compute node to node bed slope using forward, backward and central differences

    for j = 1:Nn
    
        if j == 1

            S(1,j) = ((n1(1,j) - n1(1,j+1)) / dx);

        elseif j > 1 && j <= Nn - 1

            S(1,j) = ((n1(1,j-1) - n1(1,j+1)) / (2*dx));          

        else

            S(1,j) = ((n1(1,j-1) - n1(1,j)) / dx);

        end

        % Compute water depth using manning-strickler resistance formula

        if S(1,j) > 0

            y(1,j) = ((Ksx(1,j)^0.33*qw(1,j)^2)/(a^2*g*S(1,j)))^0.30;

        else

            y(1,j) = y(1,j-1);

        end
        
        % Compute water surface elevation
        
        WSE(1,j) = n1(1,j) + y(1,j);
        
        % Compute the velocity based on the contiuity equation

        V(1,j) = qw(1,j) / y(1,j);

        % Compute dimensional shear stress distribution per node to node slope

        Tbed(1,j) = Rho * g * y(1,j) * S(1,j);

        % Compute dimensionless shear stress distribution per node to node slope - only works for steady flow

        Tbede(1,j) = ((Ksx(1,j)^0.33 * qw(1,j)^2) / (a^2 * g))^0.30 * ((S(1,j))^0.70 / (Rd * Dgx(1,j)));

        % Compute the shear velocity

        Ushear(1,j) = ((Ksx(1,j)^0.33 * qw(1,j)^2) / (a^2))^0.15 * S(1,j)^0.35 * g^0.35;
    
    end

end

