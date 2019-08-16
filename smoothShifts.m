function [dSA] = smoothShifts(dS, LD, k)
    %Smooth all shifts: average of k-shifts and k- points
    dSA = dS;
    for i = 2:length(dS) - 1
        p = 0; ds = 0;
        n = 0;
        
        %Use a point with non-zero shift
        if (norm(dS(i, :)) > 0)
            for j = i - k : i + k
                %Jump points with zero shifts
                if (j > 0) && (j <= length(dS)) && (norm(dS(j, :)) > 0)
                        p = p + LD(j, :);
                        ds = ds + dS(j, :);
                        n = n + 1;
                end
            end

            %Average coordinate of (LD(i-k), LD(i+k)) and average shift of (dS(i-k), dS(i+k))
            p = p / n;
            ds = ds / n;

            %New shift: sum of both
            dSA(i, :) = p - LD(i, :) + ds;
        end
    end
end

