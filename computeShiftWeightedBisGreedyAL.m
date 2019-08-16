function [ds, dmin, res] = computeShiftWeightedBisAL(AL, LD, Ls, NP, idx1, i, dmax, pos_bar_disp, pos_disp_bar, pos_bar_disp_next)
    %Interpolate direction vectors sn for all intermediate points of the segment
    dmin = dmax;
    ds = 0;
    res = false;
    
    %Get distance of LD(i) to the nearest point of L
    lmin = NP(i, 1); imin = NP(i, 2); dmin = NP(i, 3); pmin = NP(i, 4:5);
    
    %Shift vertex by the ratio
    if (dmin < dmax)

        %Point inside outer buffer, compute r
        r = 1 - dmin ./ dmax;

        %Get vector LN(i) - pmin
        u = LD(i, :) - pmin;
        un = u / sqrt(sum(u'.^2,1));

        %Compute interpolated shift directions
        s = Ls(min(i - idx1 + 1, length(Ls)), :) - LD(i, :);
        ns = sqrt(sum(s'.^2,1));
        sn = s / ns;
        %plot([LD(i, 1); LD(i, 1) + s(1, 1)], [LD(i, 2); LD(i, 2) + s(1, 2)], 'r');

        %Compute angle A2 between un and sn
        arg = dot(un, sn) / (norm(un) * norm(sn));
        arg = min(arg, 1); arg = max(arg, -1);
        A2 = acos(arg);
        
        %Compute shift ratio in A2 direction
        if (abs(A2 - pi/2) < 0.0001)
            R = r;
        else
            R = r / cos(A2);
        end

        %Project point LD[i] on the buffer in A2 direction
        Q = LD(i, :) + (R + 0.001) * dmax * sn;
       
%         B = polybuffer(AL{lmin},'lines', dmax);
%         I = intersect(B, [LD(i, :); Q]);
%         plot(B, 'FaceColor', 'y', 'EdgeColor', 'y', 'FaceAlpha', dmax);
%         I(any(isnan(I), 2), :) = [];
        
        %Create part of the buffer: line segment between L[imin], L[imin+1] given by Q1, Q2
        %Shift by bisector to the right halfplane
        b = bisector(AL{lmin}(imin, :), 0.5 * (AL{lmin}(imin, :) + AL{lmin}(imin + 1, :)), AL{lmin}(imin + 1, :));
        bn = b / sqrt(sum(b'.^2,1));
        
        %Switch orientation of the normal to the left halfplane
        if (pos_disp_bar == 0)
            bn = - bn;
        end
        
        Q1 = AL{lmin}(imin, :) + dmax * bn;
        Q2 = AL{lmin}(imin + 1, :) + dmax * bn;
        %plot([Q1(1, 1); Q2(1, 1)], [Q1(1, 2); Q2(1, 2)], 'k');

        %Test intersection with:
        %   cirle at L[i]
        %   segment L[i], L[i+1]
        %   cirle at L[i+1]
        I1 = []; I2 = []; I3 = [];

        I1 = getLineCircleIntersection(LD(i, 1), LD(i, 2), Q(1,1), Q(1, 2), AL{lmin}(imin, 1), AL{lmin}(imin, 2), dmax);
        I2 = get2LineSegmentsIntersection(LD(i, 1), LD(i, 2), Q(1,1), Q(1, 2), Q1(1, 1), Q1(1, 2), Q2(1, 1), Q2(1, 2));
        I3 = getLineCircleIntersection(LD(i, 1), LD(i, 2), Q(1,1), Q(1, 2), AL{lmin}(imin + 1, 1), AL{lmin}(imin + 1, 2), dmax);
        I = [I1; I2; I3];

        %Shift position
        pos_shift = 1 - pos_bar_disp;
        
        %Find the most suitable intersection point on the buffer from I1 - I3
        %It has the same position according to (LD(i), LD(i+1))
        [mi, ni] = size(I);
        rmin = 1; kmin = -1;
        for k = 1 : mi
            %Compute r(I)
            [imink, dmink, pmink, cmink] = findNearestPoint2(I(k, :), AL{lmin});
            ri = abs(1 - dmink ./ dmax);

            %Position of I according to (LD(i), LD(i+1))
            if (i == length(LD))
                pos_int_disp = getPointLinePosition(I(k, :), LD(i - 1, :), LD(i, :));
            else
                pos_int_disp = getPointLinePosition(I(k, :), LD(i, :), LD(i + 1, :));
            end

            %Intersection I on the same side of (LD(i), LD(i+1)) as the shift position 
            if (pos_int_disp == pos_shift)
                if (ri < rmin)
                    rmin = ri;
                    kmin = k;
                end
            end
        end

        %We found point on the outer buffer in the direction of sn
        nS = 2 * dmax;
        if (rmin < 1.0e-3)
            %Compute the shift and ratio R
            S = I(kmin, :) - LD(i, :);
            nS = sqrt(sum(S'.^2,1));
            Sn = S / nS;
            R = nS / dmax;

            %Shift to the outer buffer
            ds = R * dmax * Sn;
            
            res = true;
            
        %Otherwise, use the direction to the shortest point of L
        else
            ds = r *dmax * sn;
            res = true;
        end
        
        %The nearest barrier between 2 point has been changed
        %Use the closer solution: interpolated point between (ps, pe) vs point on the buffer
        if (pos_bar_disp ~= pos_bar_disp_next)
            if ns < nS
                ds = s;
                
                res = true;
            end
        end
    end
end

