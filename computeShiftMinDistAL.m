function [ds, dmin] = computeShiftMinDistAL(LD, NP, i, dmax)
    %Interpolate direction vectors sn for all intermediate points of the segment
    dmin = dmax;
        
    %Get distance of LD(i) to the nearest point of L
    %[imin, dmin, pmin, cmin] = findNearestPoint2(LD(i, :), L);
    imin = NP(i, 2); dmin = NP(i, 3); pmin = NP(i, 4:5);
        
    %Shift vertex by the ratio, if necessary
    if (dmin < dmax)
        
        %Compute ratio
        r = 1 - dmin ./ dmax;

        %Compute shift vector
        u = (LD(i, :) - pmin);
        un = u / sqrt(sum(u'.^2,1));
                
        %Shift to the outer buffer
        ds = r * dmax * un;
    end
end

