function [dS] = applyPropagation(dS, ten)
   %Apply propagation to first and last points
    K = ones(length(dS), 1);
    
    %Compute the propagation coefficient
    for i = 1 : length(dS)
        kp = propagationCoeff(i, 1, length(dS), ten);
        K(i) = kp;
    end
    
    %Apply the propagation
    dS = K .* dS;
end

