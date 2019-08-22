function [imin, dmin, pmin, cmin, points, dist] = findPointsCloserThan(p, L, dth)
        %Find point on the polyline L nearest to p

        %Trim L: remove last point
        LT = L(1:length(L)-1,:);

        %Project Pi to Pj, Pj+1
        u1 = diff([L(:, 1), L(:, 2)]);
        u2 = p - LT;

        %Dot product and the norm of the rows
        dij = dot(u1', u2');
        nij = sqrt(sum(u1'.^2,1));
        c = dij./(nij .* nij);

        %Intersection outside the line segment
        %Move it to the closest vertex
        c(c>1) = 1; c(c<0) = 0;

        %Create intersection point Pp
        Pp = [LT(:, 1) LT(:,2)] + c' .* u1;

        %Distance d(i, Pp)
        u3 = p - Pp;
        d = sqrt(sum(u3'.^2,1));

        %Find point with minimum distance
        [dmin, imin]= min(d);
        pmin = Pp(imin,:);

        %Find indices of points closer than dth
        idxs = find(d < dth);
        points = Pp(idxs, :);
        dist = (d(idxs))';
        
        %Closest point: parameter c (internal point or vertex?)
        cmin = c(imin);
end

