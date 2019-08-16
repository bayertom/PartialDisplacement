function [posit] = getPointPolylinePosition(p, L)
      %Get distance of p to the closest point of L
      [imin, dimin, pimin, cimin] = findNearestPoint2(p, L);
     
      %Get its position
      posit = getPointLinePosition(p, L(imin, :), L(imin + 1, :));
end

