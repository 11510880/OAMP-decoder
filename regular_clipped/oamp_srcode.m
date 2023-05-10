clear all;clc;
% constants
B = 64;
L = 2048;
CR = -13;
SNR = 2;
R = 0.5;
maxIter = 50;
N = B * L;

% generate sparse x
x = generate_sparse_verctor(B,L);
M = L * log2(B) / R;
z_full = dct(x);
rng("default");
order = randperm(N, M);
z = z_full(order);
% get clipped threshold
epsilon = get_threshold_by_cr(CR,z);
signal = clip(z, epsilon);
[y, alpha, sigma] = linear_model_clip(signal,SNR);   

% initialize the input to the LMMSE part
zPri = zeros(N, 1); 
vZpri = 1;
delta = M / N;
% output variance of estimation of NLE, i.e., prediction of MSE
vzpris = zeros(1, maxIter);
% output variance of estimation of LE
vxpris = zeros(1, maxIter);
% the true mse each iteration
mses = zeros(1, maxIter);
for i=1:maxIter
   fprintf("iteration number %d \n", i);
   vzpris(i) = vZpri;
   [zPost,vZpost] = declip(y, zPri(order), vZpri, sigma, epsilon);
   % zPost M * 1,  zPri, N * 1
   [vZorth, zOrth] = orthogonalization(vZpost, vZpri, zPost, zPri, order);
   
   vXpri = vZorth;
   vxpris(i) = vXpri;
   xPri = idct(zOrth);
   % vXpri = max(vXpri, 1e-6);    
   [xPost,vXpost] = sr_demodulation(B, L, vXpri, xPri);
   % xPost N * 1,  xPri, N * 1
   [vXorth,xOrth] = orthogonalization(vXpost,vXpri, xPost, xPri);
   mses(i) = norm(xOrth - x)^2 / length(x);
   vZpri = vXorth;
   zPri = dct(xOrth);
   fprintf("vXpri: %d, vZpri: %d, mse: %d \n", vXpri, vZpri, mses(i));
end


figure
semilogy(1:maxIter, vzpris, "-x");
hold on
semilogy(1:maxIter, mses, '-^');
legend("Predicted MSE", "True MSE");
xlabel("iteration number");
ylabel("MSE");

figure
semilogx(vzpris, vxpris, 'r');
hold on
semilogx([vzpris(2:end), vZpri], vxpris,  'b');
legend("de-clip", "de-modulation");
xlabel("vz");
ylabel("vx");






