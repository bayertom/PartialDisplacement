function [LDI, dSI] = removeIntersectedDirections(LD, dS)
    %Remove intersected directions
    %Sort dS by length
    nDs = sqrt(sum(dS'.^2,1));
    dS(:, 3) = nDs';
    dS(:, 4) = 1:length(dS);
    %dSs = sortrows(dS, 3,'descend');
    dSs = sortrows(dS, 3);
    %dSs = dS;
    
    %Add the third column, if a point generates direction with no
    %intersections
    LD(:, 3) = ones(length(LD), 1);
    
    %Assign
    %dSs = dS;

    %Process all segments
    for i = 1 : length(dS)
        
        %Analyze non-zero shift vectors
        if (dSs(i, 3) > 0)
            
            %Get index of the segment
            idx = dSs(i, 4);
            
            %Create line segment
            pi = LD(idx, 1:2);
            pii = pi + dS(idx, 1:2);

            %Test acoording the remaining segments
            intersect = false;
            for j = 1 : length(dS)
                if (j ~= idx) && (LD(j, 3) == 1)
                    
                    %Create line segment
                    pj = LD(j, 1:2);
                    pjj = pj + dS(j, 1:2);

                    %Find intersection of both segments
                    I = get2LineSegmentsIntersection(pi(1, 1), pi(1, 2), pii(1, 1), pii(1, 2), pj(1, 1), pj(1, 2), pjj(1, 1), pjj(1, 2));

                    %Intersection has been found
                    if (length(I) > 0)
                        LD(idx, 3) = 0;
                        break;
                    end
                end
            end
            
%             %Direction Sn(i) has no intersection, mark point as usable
%             if (intersect == false)
%                 LD(idx, 3) = 1;
%             end
        end
    end
    
    %Copy all points generating Sn without intersections
    LDI = []; dSI = [];
    for i = 1 : length(LD)
        if (LD(i, 3) == 1)
            LDI = [LDI; LD(i, 1:2)];
            dSI = [dSI; dS(i, 1:2)];
        end
    end
end

