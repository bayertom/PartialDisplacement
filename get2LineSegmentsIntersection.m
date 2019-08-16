function [I] = get2LineSegmentsIntersection(x1, y1, x2, y2, x3, y3, x4, y4)
%Get intersection of 2 lines
I = [];

%Compute denominator
denom =	(y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1) ;

%Segments are parallel
if ( abs (denom) < 0.0001 )
    return;
end

%Compute numerators
numer1 = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3));
numer2 = ((x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3));

%Compute parameters s,t
s = numer1 / denom;
t = numer2 / denom;

%Segments do not intersect: do not compute any intersection
if ((s < 0.0) || (s > 1) || (t < 0.0) || (t > 1))
    return;
end

%Compute intersection
xi = x1 + s * (x2 - x1);
yi = y1 + s * (y2 - y1);

I = [xi, yi];

end

