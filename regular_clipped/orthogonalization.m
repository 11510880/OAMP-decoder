function [vOrth,uOrth] = orthogonalization(vPost,vPri, uPost, uPri, order)
%ORTHOGONALIZATION  normalization of means and variances
if nargin == 5
    M = length(uPost);
    N = length(uPri);
    delta = M /N;
    vPost_ = delta * vPost + (1 - delta) * vPri;
    vOrth = 1 / (1./vPost_ - 1./vPri);
    vOrth = max(1e-8, vOrth);
    uPost_ = zeros(N,1);
    uPost_(order) = uPost;
    uPost_(setdiff(1:N,order)) = uPri(setdiff(1:N,order)); 
    uOrth = vOrth .* (1./vPost_ .* uPost_ - 1 / vPri .* uPri);
else
    vOrth = 1 / (1./vPost - 1./vPri);
    uOrth = vOrth .* (1./vPost .* uPost - 1 / vPri .* uPri);
end
end

