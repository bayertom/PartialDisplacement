function [LN, IDXN] = densifyPolyline(L, ds)
    %Densify polyline segments with a given length step ds
    %New points are equally spaced
    LN =[];
    IDXN = [];

    %Compute differences
    u = diff([L(:, 1), L(:, 2)]);

    %Normalize vector
    nu = sqrt(sum(u'.^2,1));
    un = (u' ./ nu)';

    %Densify polyline segments one by one 
    eps = 0.0001;
    for i = 1:length(L) - 1

        %Compute parameter alpha, direction and densify
        if (ds > 0)

            %Index of old points in the densified polyline
            IDXN = [IDXN; 1 + length(LN)];

            alpha = 0:ds:nu(i) - eps;
            LNi = L(i, :) + alpha' * un(i, :);

            %Merge to the previous segment
            LN = [LN; LNi];
            
        %No densification
        else
            IDXN = [IDXN; i];
            LN = [LN; L(i, :)];
        end
        
         %Add end point of the last segment
         if (i == length(L) - 1)
             LN = [LN; L(length(L), :)];
         end
    end
    
    %Add last vertex
    IDXN = [IDXN; length(LN)];
end



