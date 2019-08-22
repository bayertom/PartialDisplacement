function [k] = propagationCoeff(j,  idxs, idxe, ten)
    %Compute the propagation coefficient k in the direction  vector sn
    k = 1;

    %Propagate points inside the outer buffer and outside the inner buffer
    %Propagate at least 1/3 of the  first and last points
    nprop = max(floor(abs(ten) * (idxe - idxs) / 3), 3);

    %At least 2 points will be propagated
    if (nprop > 2)
        %Displace n initial points: k(idxs) = 0, k(idxs + n) = 1
        if (j < nprop + idxs )
            k = 1 + ten * ((idxs + nprop - j) / nprop)^1;

        %Displace n last points: k(idxe) = 0, k(idxe - n) = 1
        elseif (j > idxe - nprop)
            k = 1 + ten * ((j-idxe + nprop ) / nprop)^1;
        end
    end
end

