function [idxl, idxmin, dmin, pmin, cmin] = findNearestPointAL(p, AL)
    %Find point on the list of polylines AL nearest to p
    idxl = - 1; idxmin = -1; dmin = 9999; pmin = []; cmin = [];
    
    %Process all lines in the list
    for i = 1 :length(AL)
        [imin, dimin, pimin, cimin] = findNearestPoint2(p, AL{i});
        
        if (dimin < dmin)
            idxl = i;           %Index of the line
            idxmin = imin;      %index of the point
            dmin = dimin;       %Nearest distance
            pmin = pimin;       %Nearest point
            cmin = cimin;       %Parameter c: vertex of L or intermediate point?
        end
    end
end
