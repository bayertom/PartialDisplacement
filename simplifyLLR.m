function [LS] = simplifyLLR(L, llr_min, perc)
%Simplify line using the LLR criterion
LLR = ones(length(L), 2);
IDX = zeros(length(L), 1);

%Compute LLR for all points
for i = 2 : length(L) - 1
    llr = getLLR(L(i-1, 1), L(i-1, 2), L(i, 1), L(i, 2), L(i+1, 1), L(i+1, 2));
    LLR(i, :) = [i, llr];
end

%Sort LLR in descendent order
LLR = sortrows(LLR, 2, 'descend');

%Mark all important vertices
i = 1;
while (LLR(i, 2) > llr_min)
     IDX(LLR(i, 1)) = 1; 
     i = i + 1;
end

%Keept at least percentage of less important points
n_dest = min(perc, 0.5) * (length(L) - i);

%If all important points found, use some less important
is = i;
while (i < length(L)) && (LLR(i, 2) < llr_min) && (i - is < n_dest)
    IDX(LLR(i, 1)) = 1; 
    i = i + 1;
end

%Create simplified line
LS = [];
LS = [LS; L(1,:)];

%Add all important vertices
for i = 2 : length(L) - 1
    if IDX(i) == 1
        LS = [LS; L(i, :)];
    end
end

%Add last point
LS = [LS; L(length(L), :)];