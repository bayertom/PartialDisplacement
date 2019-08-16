function [dir] = findShiftDirectionAL(AL, LD, IDX)
    %Find position of pmin against LD (barrier against displaced)
    %Return first value different from -1
    dir = 0;
    
    %Process all barriers
    for i = 1 :length(AL)
    
        %Process all original vertices of LD before the densification
        %Compute the initial direction vector sn
        nj = length(IDX);
        for j = 1 : nj

            %Get distance of LD(j) to the closest point of L
            [jmin, djmin, pjmin, cjmin] = findNearestPoint2(LD(IDX(j), :), AL{i});

            %Position of pmin against LD; the closest point is known (barrier against displaced)
            pos_bar_disp = getPointPolylinePosition(pjmin, LD(IDX, :));

            %First shift direction has been fond (left or right)
            if (pos_bar_disp ~= -1)
                dir = pos_bar_disp;
                return;
            end
        end
    end
end

