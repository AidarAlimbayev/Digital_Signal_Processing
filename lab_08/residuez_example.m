% Example:
    %   Compute the partial fraction expansion of the following transfer
    %   function H(z) = (1 + 2z^-1) / (1 - z^-1 + 2z^-2).
 
    num = [1 1];                % Numerator coefficients
    den = [1 -1 2];             % Denominator coefficients
    [r,p] = residuez(num,den)   % H(z) = r(1)/(1-p(1)z^-1) + ...
                                %        r(2)/(1-p(2)z^-1)