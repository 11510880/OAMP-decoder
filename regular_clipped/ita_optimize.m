function [resultDenoised] = ita_optimize(tao, r, sparseP)
% here mmseB is an average value computed outside, not a function
% ita_mmse is a function
% tao is the estimation of noise, i.e. r = x + N(0, tao)
% sparseP is the ratio of non-zero entries of original signal x

mPri = 0;
vPri = 1/sparseP;
% (mPri, vPri, sparseP) is the parameter of Gaussian part for the spike and
% slab prior (Bernoulli Gaussian)
[uPost,vPost] = ita_mmse(tao, r, vPri, mPri, sparseP);
mmseB = mmse_b(vPost);

C = tao / (tao - mmseB);
resultDenoised = C .* (uPost - mmseB / tao .* r);
%resultDenoised = uPost;
end

