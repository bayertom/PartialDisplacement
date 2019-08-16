function [idxl, idxmin, dmin, pmin, cmin, al_points, al_dist] = findPointsCloserThanAL(p, AL, dth)
    %Find point on the list of polylines AL nearest to p
    idxl = - 1; idxmin = -1; dmin = 9999; pmin = []; cmin = []; al_points = []; al_dist = [];
    
    %Process all lines in the list
    for i = 1 :length(AL)
        [imin, dimin, pimin, cimin, points, dist] = findPointsCloserThan(p, AL{i}, dth);
        al_points = [al_points; points];
        al_dist = [al_dist; dist];
        
        if (dimin < dmin)
            idxl = i;           %Index of the line
            idxmin = imin;      %index of the point
            dmin = dimin;       %Nearest distance
            pmin = pimin;       %Nearest point
            cmin = cimin;       %Parameter c: vertex of L or intermediate point?
        end
    end
end
