function [bn] = bisector(p1, p2, p3)
    %Compute normalized disector in p2 oriented to the right half-plane
    u = p1 - p2;
    v = p3 - p2;
    
    %Norms of vectors
    nu = sqrt(sum(u'.^2,1));
    nv = sqrt(sum(v'.^2,1));
    
    %Bisector
    b1 = u / nu + v / nv;
    b = nu * v + nv * u;
    nb = sqrt(sum(b'.^2,1));
    
    %Colinear vectors
    if (nb < 1.0e-6)
        b = [-u(:, 2), u(:, 1)] / nu;
    end
    
    %Switch orientation
    T = [u; b];
    t = det(T);
    if (t < 0)
        b = -b;
    end
    
    %Normalize bisector
    nb = sqrt(sum(b'.^2,1));
    bn = b / nb;
   
end

