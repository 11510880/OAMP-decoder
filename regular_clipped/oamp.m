% main entrance of Orthogonal AMP
clear all
clc
N = 4000;
M = 2000;
sparseP = 0.2;
kappa = 5;
SNR = 60;
maxIter = 50;

% generate source signal
signal = bernoulli_gaussian_generator(sparseP, N)';
[A, Lambda] = ill_condition_matrix_generator(kappa, N, M);
measurement = A * signal;
sigma = norm(measurement) .^2 / M / 10 ^ (SNR / 10);
noise = sqrt(sigma) * randn(M, 1);
% fprintf("snr %d", snr(measurement, noise));
y = measurement + noise;

% plot SE
v_t = 1 / N * (norm(signal) .^2);
se_v = zeros(1, maxIter);
se_tao = zeros(1, maxIter);
for i = 1:maxIter
    se_v(i) = v_t;
    mmseA = mmse_a(sigma, v_t, Lambda);
    tao_t = se_le(v_t, mmseA);
    se_tao(i) = tao_t;
    r = signal + sqrt(tao_t) * randn(N, 1);
    [uPost,vPost] = ita_mmse(tao_t, r, 1/sparseP, 0, sparseP);
    mmseB = mmse_b(vPost);
    v_t = se_nle(mmseB, tao_t);
    v_t = max(v_t, 0.0001);
end
figure
semilogy(1:maxIter, se_v, '-^');
hold on 
semilogy(1:maxIter, se_tao, '-*');
legend("se v", "se tao ");


% initialization
s = zeros(N, 1);
v_hats = zeros(1, maxIter);
tao_hats = zeros(1, maxIter);
mse_s = zeros(1, maxIter);
v_t = 0;
for iter = 1:maxIter
    %fprintf("iteration %d mse before mmse: %d\n", iter, 1/ N * (norm(s - signal) .^2));
    v_hat = max(0.0000001, (norm(y - A * s) .^2 - M * sigma) / trace(A'* A));
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
    % predict using state evolution
    if iter == 1
        v_t = v_hat;
    end
    se_v(iter) = v_t;
    mmseA = mmse_a(sigma, v_t, Lambda);
    tao_t = se_le(v_t, mmseA);
    se_tao(iter) = tao_t;
    [uPost,vPost] = ita_mmse(tao_t, r, 1/sparseP, 0, sparseP);
    mmseB = mmse_b(vPost);
    v_t = se_nle(mmseB, tao_t)     ;
    v_t = max(v_t, 0.0000001);
    fprintf("tao_se: %d, v_se:%d \n", tao_t, v_t);
end
figure
semilogy(1:maxIter, mse_s, '-*'); hold on;
semilogy(1:maxIter, v_hats, '-^');           
semilogy(1:maxIter, se_v, '-.');
legend("Simulated MSE", "Etismated MSE", "Predict MSE");

% figure
% plot(1:maxIter, mse_s, '--');
% hold on
% plot(1:maxIter, v_hats, '-^');
% % hold on
% % plot(1:maxIter, v_t, '-*');
% legend("Simulated MSE", "Estimated v_t")
% title("Simlated MSE and Prediction");



