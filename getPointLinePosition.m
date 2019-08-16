function [posit] = getPoinLinePosition(p, p1, p2)
      
      %Get position of the point a and line l(p1, p2)
      t = det([p - p1; p2 - p1]);
      
      %Point in the right halfplane
      if (t > 0)
          posit = 1;
      
      %Point in the left halfplane
      elseif (t < 0)
          posit = 0;
          
      %Point on the line
      else
          posit = -1;
      end
end

