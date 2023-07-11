% main entrance of Orthogonal AMP
clear all
clc
N = 1500;
M = 2000;
sparseP = 0.2;
kappa = 5;
SNR = 60;
maxIter = 35;

% generate source signal
signal = bernoulli_gaussian_generator(sparseP, 3 * N)';
A_original = randn(M, N);

A = zeros(M, 3 * N);
% build lattice
for i=1:N
    A(:,(i-1) * 3 +1)= A_original(:,i);
    A(:,(i-1) * 3 +2)= 2 * A_original(:,i);
    A(:,(i-1) * 3 +3)= 3 * A_original(:,i);
end

N = 3 * N;
% [A, Lambda] = ill_condition_matrix_generator(kappa, N, M);
measurement = A * signal;
sigma = norm(measurement) .^2 / M / 10 ^ (SNR / 10);
noise = sqrt(sigma) * randn(M, 1);
% fprintf("snr %d", snr(measurement, noise));
y = measurement + noise;

% initialization
s = zeros(N, 1);
v_hats = zeros(1, maxIter);
tao_hats = zeros(1, maxIter);
mse_s = zeros(1, maxIter);
v_t = 0;
for iter = 1:maxIter
    %fprintf("iteration %d mse before mmse: %d\n", iter, 1/ N * (norm(s - signal) .^2));
    v_hat = max(1e-6, (norm(y - A * s) .^2 - M * sigma) / trace(A'* A));
    wOptimize = w_optimize(A, v_hat, sigma);
    r = s + wOptimize * (y - A * s);
    B = eye(N) - wOptimize * A;
    tao_hat = 1 / N * (trace(B * B') * v_hat + trace(wOptimize * wOptimize') * sigma);
    s = ita_optimize(tao_hat, r, sparseP);
    mse = 1 /N * (norm(s - signal) .^2);
    mse_s(iter) = mse;
    tao_hats(iter) = tao_hat;
    v_hats(iter) = v_hat;
    fprintf("iteration %d \n", iter);
    fprintf("h = r - x, E{h^2} is  %d, estimated tao is: %d \n", 1 /N * (norm(r - signal) .^2), tao_hat)
    fprintf("q = s - x, E{q^2} is: %d, estimated v is %d \n", mse, v_hat);
end
figure
semilogy(1:maxIter, mse_s, '-*'); hold on;
semilogy(1:maxIter, v_hats, '-^');           
legend("Simulated MSE", "Etismated MSE");



