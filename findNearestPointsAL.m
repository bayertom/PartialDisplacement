function [idxmin, dmin, pmin, cmin] = findNearestPointsAL(p, AL)
    %Find the nearest point of the each polyline on the list to p
    idxmin = []; dmin = []; pmin = []; cmin = [];
    
    %Process all lines in the list
    for i = 1 :length(AL)
        [imin, dimin, pimin, cimin] = findNearestPoint2(p, AL{i});
        
        %Add to the list
        idxmin = [idxmin; imin];
        dmin = [dmin; dimin];
        pmin = [pmin; pimin];
        cmin = [cmin; cimin];
     
    end
end
