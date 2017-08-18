function [ Qwseries ] = Hydrograph( dt,T,tt,Qw,HYDTime,Qwincremental )

% Function to compute a hydrograph array from data at select points in time
  
Qwseries = zeros(tt,1);
Qwvolume = zeros(tt,1);
Qwavgseries = zeros(tt,1);
        
    for i = 1:tt

        if i == 1

            Qw1(1) = Qw;

        elseif i > 1 && i <= HYDTime(2,1)

            Qw1(1) = Qw1 + Qwincremental(1,1);

        elseif i > HYDTime(2,1) && i <= HYDTime(3,1)

            Qw1(1) = Qw1 + Qwincremental(2,1);

        elseif i > HYDTime(3,1) && i <= HYDTime(4,1)

            Qw1(1) = Qw1 + Qwincremental(3);

        elseif i > HYDTime(4,1) && i <= HYDTime(5,1)

            Qw1(1) = Qw1 + Qwincremental(4,1);

        elseif i > HYDTime(5,1) && i <= HYDTime(6,1)

            Qw1(1) = Qw1 + Qwincremental(5,1);

        elseif i > HYDTime(6,1) && i <= HYDTime(7,1)

            Qw1(1) = Qw1 + Qwincremental(6,1);

        elseif i > HYDTime(7,1) && i <= HYDTime(8,1)

            Qw1(1) = Qw1 + Qwincremental(7,1);

        elseif i > HYDTime(8,1) && i <= HYDTime(9,1)

            Qw1(1) = Qw1 + Qwincremental(8,1);

        elseif i > HYDTime(9,1) && i <= HYDTime(10,1)

            Qw1(1) = Qw1 + Qwincremental(9,1);  
        
        elseif i > HYDTime(10,1) && i <= HYDTime(11,1)

            Qw1(1) = Qw1 + Qwincremental(10,1);      
        
        elseif i > HYDTime(11,1) && i <= HYDTime(12,1)

            Qw1(1) = Qw1 + Qwincremental(11,1);       
            
        end

    Qwseries(i) = Qw1(1);

    Qwvolume(i) = Qw1(1) * dt;
    
    end
    
    Qwsum = sum(Qwvolume);
    
    Qwavg = Qwsum / T;
    
    for i = 1:tt;
    
        Qwavgseries(i) = Qwavg;
    
    end
    
    save('Qwseries.mat','Qwseries');
    save('Qwaverage.mat','Qwavg');
    save('Qwaverageseries.mat','Qwavgseries');
        
end

