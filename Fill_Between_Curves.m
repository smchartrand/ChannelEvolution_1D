function [ output ] = Fill_Between_Curves(StationCrop,DEMProfileAdjust,ZeroCrossFunc)
% This function finds the exact intersection points between two curves in
% order to use the fill function to fill between two curves.
% The function was written by Mike Garrity of MATLAB - 
% http://blogs.mathworks.com/graphics/2015/10/13/fill-between/
    x = StationCrop;
    y1 = DEMProfileAdjust;
    y2 = ZeroCrossFunc;        
    output = [];
    % Calculate t in range [0 1]
    calct = @(n) (n(3,1)-n(2,1))/(n(3,1)-n(2,1)-n(3,2)+n(2,2));
    % Generate interpolated X and Y values
    interp = @(t,n) n(:,1) + t*(n(:,2) - n(:,1));
    for i=1:length(x)
    % If y2 is below y1, then we don't need to add this point.
        if y2(i) <= y1(i)
            % But first, if that wasn't true for the previous point, then add the
            % crossing.
            if i>1 && y2(i-1) > y1(i-1)
                neighborhood = [x(i-1), x(i); y1(i-1), y1(i); y2(i-1), y2(i)];
                t = calct(neighborhood);
                output(:,end+1) = interp(t,neighborhood);
            end
        else
        % Otherwise y2 is above y1, and we do need to add this point. But first
        % ...
            % ... if that wasn't true for the previous point, then add the
            % crossing.
            if i>1 && y2(i-1) <= y1(i-1)
                neighborhood = [x(i-1), x(i); y1(i-1), y1(i); y2(i-1), y2(i)];
                t = calct(neighborhood);
                output(:,end+1) = interp(t,neighborhood);
            end

            % add this point.
            output(:,end+1) = [x(i); y2(i); y1(i)];
        end
        
    end

end

