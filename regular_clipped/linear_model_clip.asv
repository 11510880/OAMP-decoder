function [y, alpha, sigma_] = linear_model_clip(signal, SNR)
%LINEAR_MODEL_CLIP linear model with a clipped verison of dct input
% y = alpha * clip(signal,epsilon) + noise, where signal = Fx and x is a sparse vector, F is
% dct-matrix, noise ~ N(0, sigma * I)
alpha = 1 / sqrt(norm(signal) .^2 / length(signal));
% signalNormalized = signal .* alpha;

sigma = 1 / 10 ^ (SNR / 10);
% sigma = norm(signalNormalized) .^2 / length(signal) / 10 ^ (SNR / 10);
noise = sqrt(sigma) / alpha * randn(length(signal), 1);
% y = signalNormalized + noise;
y = signal + noise;
% y = y / alpha;
% sigma_ = sigma;
sigma_ = sigma / (alpha.^2);
end

