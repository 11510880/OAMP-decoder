function [wOptimize] = w_optimize(A, v, sigma, W)
% a LMMSE optimized linear MMSE Estimator
% OMAP 15c, notice that all v, tau,... represents v * v, tau * tau in the
% paper
if nargin == 3
   w_lmmse = @(v,sigma,A)v * A'/(v * A * A' + sigma * eye(size(A, 1)));
   W = w_lmmse(v, sigma, A); 
end
% OMAP lemma 1, 38a
N = size(A,2);
wOptimize = N / trace(W * A) * W; 
end
