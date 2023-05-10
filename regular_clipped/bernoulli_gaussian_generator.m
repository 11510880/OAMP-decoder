function [sparseArray] = bernoulli_gaussian_generator(p, N)
%BERNOULLIGAUSSIAN_GENERATOR generate Bernoulli-Gaussian variables, the
%variance of Gaussian is 1/p

sparseArray = binornd(1, p, [1,N]);
for i = 1:N
    if sparseArray(i) ~= 0
        sparseArray(i) = sqrt(1/p) * randn(1);
    end
end
end

