function [uPost,vPost] = ita_mmse(tao, r, vPri, mPri, sparseP)
%ITA_MMSE MMSE denoiser
% R = X + N(0, tao), X is a bernoulli Gaussian variable with p as sparsity
% the input is r, a vector form of R, all operations are elementwise
% level (1-p sparse postions)
% return uPost:E(X|R), vPost:Var(E|R)
p = sparseP;
ug = mPri;
vg = vPri;
v = tao;
y = r;
alpha = sqrt((v + vg) / v);
beta = 1/2 * ((y - ug).^2 / (v + vg) - y.^2 / v);
pPrime = p./ (p + (1-p) .* alpha .* exp(beta));
vgPrime = 1 / (1 / vg + 1 / v);
ugPrime = (1 / vg * ug + 1 / v * y) * vgPrime;

uPost = pPrime.* ugPrime;
vPost = (pPrime - pPrime.^2) .* ugPrime + pPrime .* vgPrime;
end

