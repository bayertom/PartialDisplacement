function [llr] = getLLR(x1, y1, x2, y2, x3, y3)
    %Compute Local Length Ratio
	ux = x1 - x2;
	uy = y1 - y2;
	vx = x3 - x2;
	vy = y3 - y2;

	%Lengths of vectors
	d1 = sqrt(ux * ux + uy * uy);
	d2 = sqrt(vx * vx + vy * vy);

	%Normalize vectors
	unx = ux / d1;
	uny = uy / d1;
	vnx = vx / d2;
	vny = vy / d2;

	%Compute points A=[xa, ya], B=[xb, yb] on segments
	xa = x2 + unx;
	ya = y2 + uny;
	xb = x2 + vnx;
	yb = y2 + vny;

	%Distance between A, B
	dx = (xb - xa);
	dy = (yb - ya);

	%Compute LLR = (dist(p1, p2) + dist(p2, p3)) / dist (A, B)
	llr = 2.0 / sqrt(dx * dx + dy * dy);

end

