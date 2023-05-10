function [A, Lambda] = ill_condition_matrix_generator(kappa, N, M)
%ILL_CONDITION_MATRIX_GENERATOR generate ill condition matrix 
% kappa: condition number
% output matrix size M * N
[Q,R] = qr(randn(M));
U = Q*diag(sign(diag(R)));
[Q2,R2] = qr(randn(N));
V = Q2*diag(sign(diag(R2)));

tmp = 1;
for i = 1:M-1
    tmp = tmp + kappa .^ (i/M);
end
lambdaM = N / (tmp);
eigenvalues = zeros(1, M);
eigenvalues(M) = lambdaM;
for i = M-1:-1:1
    eigenvalues(i) = kappa .^ (1/M) * eigenvalues(i+1);
end
Sigma = diag(eigenvalues);
Sigma = padarray(Sigma,[0,N-M], 0, 'post');
A = U * Sigma * V';
Lambda = [eigenvalues, zeros(1,N-M)];
end

