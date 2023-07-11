function [xPost,vXpost] = sr_demodulation(B, L, vPri, xPri)
%SR_DEMODULATION demodulation for sr codes
% xPri is a segment with B entries of x
% both outputs are vectors
p =  zeros(B,L);
Bi = diag(sqrt(B) * ones(B,1));
xPri = reshape(xPri, B, L);
for i = 1:B
    p(i, :) = exp(-vecnorm(xPri - Bi(i, :)') .^2 ./ (2 * vPri));
end

pSum = sum(p);
p = p ./ pSum;

xPost = sqrt(B) .* p;
secondMoment = B .* p;
% vXpost = max(0, mean(mean(secondMoment - xPost.*xPost)));
% vXpost = max(1e-8, mean(mean(secondMoment - xPost.*xPost)));
vXpost = mean(mean(secondMoment - xPost.*xPost));
xPost = reshape(xPost, B * L, 1);
end


