function [mmseA] = mmse_a(sigma, vt, Lambda)
% MMSEA
% sigma: noise variance, vt, noise predicted by se_nle, Lambda: eigenvalues of A 
res = 0;
 for i = 1:length(Lambda)
    res = res + sigma * vt / (sigma + vt * (Lambda(i).^2));
end
mmseA = res / length(Lambda);
end

