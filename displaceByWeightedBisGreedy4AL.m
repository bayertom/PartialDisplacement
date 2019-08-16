function [LN] = displaceByWeightedBisGreedy4AL(AL, LD, IDX, dmax, dbuff, ten, smooth)
    %Initialize new polyline: a copy of the displaced one
    LN = LD;

    %Matrix of values of the original line: indices of vertices, sn directions
    O = [];
    
    %Find initial shift direction
    pos_bar_disp_init = findShiftDirectionAL(AL, LD, IDX);
    
    %Process all original vertices of LD before the densification
    %Compute direction vectors sn (shifting directions)
    ni = length(IDX);
    for i = 1 : ni
        %Get distance of LD(i) to the closest point of L
        %[ilmin, imin, dimin, pimin, cimin] = findNearestPointAL(LD(IDX(i), :), AL);
        [ilmin, imin, dimin, pimin, cimin, points, d] = findPointsCloserThanAL(LD(IDX(i), :), AL, dmax);

        %Position of pmin against LD (barrier against displaced)
        pos_bar_disp = getPointPolylinePosition(pimin, LD(IDX, :));
        if (pos_bar_disp == -1)
            pos_bar_disp = pos_bar_disp_init;
        end

        %Get shift position
        pos_shift = 1 - pos_bar_disp;
        
        %Position of LD(IDX) against L (displaced agains barrier)
        pos_disp_bar = getPointPolylinePosition(LD(IDX(i), :), AL{ilmin});
        
        %Compute r
        r = 1 - dimin ./ dmax;

        %Get normalized bisector bn
        if (i == 1)
            bn = bisector(LD(IDX(1), :), 0.5 * (LD(IDX(1), :) + LD(IDX(2), :)), LD(IDX(2), :));
        elseif (i == length(IDX))
            bn = bisector(LD(IDX(ni-1), :), 0.5 * (LD(IDX(ni - 1), :) + LD(IDX(ni), :)), LD(IDX(ni), :));
        else
            bn = bisector(LD(IDX(i - 1), :), LD(IDX(i), :), LD(IDX(i + 1), :));
        end
        
        %Switch orientation of bn to the left half-plane
        if (pos_shift == 0)
            bn = - bn;
        end
        
        plot([LD(IDX(i), 1); LD(IDX(i), 1) + bn(:, 1)], [LD(IDX(i), 2); LD(IDX(i), 2) + bn(:, 2)], 'b');

        %Compute un as the weighted average of vectors to the closest points
        if length(points) ~= 0
            w = (1./d')';
            ui = LD(IDX(i), :) - points();
            uin =  ui ./(sqrt(sum(ui'.^2,1)))';
            un = sum(w.*uin, 1)./sum(w,1);
            
        %No points inside the outer buffer: use the nearest point
        else

            %Get vector un = LD(i) - pmin
            u = LD(IDX(i), :) - pimin;
            un = u / sqrt(sum(u'.^2,1));
        end
        plot([LD(IDX(i), 1); LD(IDX(i), 1) + un(:, 1)], [LD(IDX(i), 2); LD(IDX(i), 2) + un(:, 2)], 'g');
       

        %Compute the shift direction sn: average of both vectors (un, t*bn)
        t = (1 - r)* dmax;
        s = 0.5 * (bn + t^(1 + (ten + 1)/2) * un);
        %s = 0.5 * (bn + un);
        sn = s / sqrt(sum(s'.^2,1));
        plot([LD(IDX(i), 1); LD(IDX(i), 1) + sn(:, 1)], [LD(IDX(i), 2); LD(IDX(i), 2) + sn(:, 2)], 'r');

        %Store values to the matrix
        O = [O; [IDX(i), sn, r, pos_bar_disp, pos_disp_bar]];
    end
    
    %Find point chains inside the buffer forming the displaced parts
    %Store their endpoints and endpoints of associated line segments
    inside = false; index = 1;
    [m, n] = size(O); 
    S = [];
    NP = zeros(length(LD), 5);
    
    for i = 1 : m - 1
        %Get indices of points of the original line
        idx1 = O(i, 1); idx2 = O(i + 1, 1);
               
        %Process all points
        for j = idx1:idx2
            
            %Get distance of LN(i) to the nearest point of L, store values
            [jlmin, jmin, djmin, pjmin, cjmin] = findNearestPointAL(LD(j, :), AL);
            NP(j, :) = [jlmin, jmin, djmin, pjmin];
            
            %Point inside the buffer
            if (djmin < dmax)
                %New start point
                if (inside == false)
                    S(index, 1) = j;
                    S(index, 3) = idx1;
                end
                
                %Point inside the buffer
                inside = true;
                
            %Point outside the buffer
            else
                %New end point
                if (inside == true)
                    S(index, 2) = j - 1;
                    S(index, 4) = idx2;
                    index = index + 1;
                end
                
                %Point outside the buffer
                inside = false;
            end
        end
    end
    
    %Last point is inside the buffer
    if (inside == true)
        S(index, 2) = j;
        S(index, 4) = O(m, 1);
    end
    
    %Create matrix of shifts
    dS = zeros(length(LD), 2); 
    dSS = dS; dSP = dS; dSA = dS;
    K = ones(length(LD), 1);
    j = 1;
    
    %Process displaced parts one by one
    [m, n] = size(S);
    idx2 = -1;
    LDO = []; dSO = [];
    
    %Copy initial undisplaced points, if exist
    if (S(1, 1) ~= 1)
        LDO = [LDO; LD(1:S(1, 1) - 1, :)];
        dSO = [dSO; zeros(S(1, 1) - 1, 2)];
    end
        
    for i = 1: m
        
        %Interval of points inside the buffer
        idx1 = S(i, 1); 
        
        %Copy all undisplaced points between displaced fragments
        if i ~= 1
            LDO = [LDO; LD(idx2 + 1:idx1 - 1, :)];
            dSO = [dSO; zeros(idx1 - 1 - idx2, 2)];
        end
        
        idx2 = S(i, 2); 
        
        %Start point of the first, end point of the last segment
        idxs = S(i, 3); 
        idxe = S(i, 4);
        
        %Find start point of the segment in O (original line L)
        while ((O(j, 1) ~= idxs) && (j < length(O)))
            j = j + 1;
        end
        
        %Initialize shift behind the inner buffer
        dsmin = dmax;
        
        %Process all displaced points inside the outter buffer
        k = idx1;
        dSk = []; LDk = [];
        while k <= idx2
            
            %We found end point of the displaced segment, increment j
            %Skip last displaced point
            if ((k == O(j + 1, 1)) && (k ~= idxe))  
                j = j + 1;
                continue;
            
            %We found start point of the segment: densify segment
            elseif ((k == O(j, 1)) || (k == idx1))              
                %Compute new points on lines of shifts
                ps = LD(O(j, 1), :) + dmax * O(j, 2:3);
                pe = LD(O(j + 1, 1), :) + dmax * O(j + 1, 2:3);
                use = pe - ps;
                dse = sqrt(sum(use'.^2,1));
                %plot([ps(:,1), pe(:, 1)],[ps(:,2), pe(:, 2)]);
                
                %Densify line segment between points of shifts
                h = dse / (O(j + 1, 1) - O(j, 1));
                [Ls, IDs] = densifyPolyline([ps; pe], h);
                
                %for hh = 1 : length(Ls)
                %    plot(Ls(hh, 1), Ls(hh, 2), 'o');
                %end
            end
            
            %Check, whether segments (LD(O,j), ps) and LD(O(j+1), pe) intersect
            I = get2LineSegmentsIntersection(LD(O(j, 1), 1), LD(O(j, 1), 2), ps(1, 1), ps(1, 2), LD(O(j + 1, 1), 1), LD(O(j + 1, 1), 2), pe(1, 1), pe(1, 2));

            if (length(I) == 0)
                
                %Compute shifts to the outer buffer for a given point LD(k)
                %Store only suitable points
                [ds, dmink, res] = computeShiftWeightedBisGreedyAL(AL, LD, Ls, NP, O(j, 1), k, dmax, O(j, 5), O(j, 6), O(j+1, 5));

                if (res)
                    dSk = [dSk; ds];
                    LDk = [LDk; LD(k, :)];
                    
                    %Find the nearest point: nearest point generates the longest shift
                    dsmin = min(dsmin, dmink);
                else
                    disp('remove')
                end
            end
            
            %Increment index of the displaced point
            k = k + 1;
        end
                
       %Apply propagation
       dSk = applyPropagation(dSk, ten);
        
       %Shift point to the inner buffer, use the longest shift generated by closest vertex
       fr = (dbuff - dsmin)/(dmax - dsmin);
       dSS = fr * dSk;
       %dSS = dSk;
       
       %Add all displaced vertices
       LDO = [LDO; LDk];
       dSO = [dSO; dSS];
         
    end
    
    %Copy final undisplaced points, if exist
    if (S(end, 2) ~= length(LD))
        LDO = [LDO; LD(S(end, 2) + 1 : length(LD), :)];
        dSO = [dSO; zeros(length(LD) - S(end, 2), 2)];
    end
   
    %Smooth shifts sn for all segments: average of k-vectors sn
    dSA = smoothShifts(dSO, LDO, smooth);
    %dSA = dSO;
    
    %Remove directions intersecting each other
    [LDI, dSA] = removeIntersectedDirections(LDO, dSA);
    
    %New position of points
    LN = LDI + dSA;
    
    %Apply LLR simplification
    LN = simplifyLLR(LN, 1.0001, 0.0);
    
    %for (i = 1 : length(LN))
    %    plot([LD(i, 1); LD(i, 1) + dSA(i, 1)], [LD(i, 2); LD(i, 2) + dSA(i, 2)], 'r');
    %end
   
end

