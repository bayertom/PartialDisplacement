function [I] = getLineCircleIntersection(x1, y1, x2, y2, xc, yc, r)
%Compute line-circle intersection
I = [];

%Reduce cordinates
x1r = x1 - xc;
y1r = y1 - yc;
x2r = x2 - xc;
y2r = y2 - yc;

%Coordinate differences and distance
dx = x2r - x1r;
dy = y2r - y1r;

%Coefficients of the quadratic equation
a = dx^2 + dy^2;
b = 2 * (dx * x1r + dy * y1r);
c = x1r^2 + y1r^2 - r^2;

%Solve quadratic equation
discr = b^2 - 4 * a * c;

%No intersection
if (discr < 0)
    return;
   
%Line is tangent, 1 intersection
elseif (discr == 0)
    
    t = -b / (2 * a);
    xi = x1 + t * dx;
    yi = y1 + t * dy;
    I = [xi, yi];

%Line is secant, 2 intersections, solve quadratic equation
else
    %First intersection
    t = (-b + sqrt(discr)) / (2 * a);
    xi1 = x1 + t * dx;
    yi1 = y1 + t * dy;
    
    %Second intersection
    t = (-b - sqrt(discr)) / (2 * a);
    xi2 = x1 + t * dx;
    yi2 = y1 + t * dy;
    
    I = [xi1, yi1; xi2, yi2];
end

