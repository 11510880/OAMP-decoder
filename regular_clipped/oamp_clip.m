function [errorSectionRate] = oamp_clip(B,L,CR,N,R,SNR,maxIter)
[x, positions] = generate_sparse_verctor(B,L);
M = L * log2(B) / R;
z_full = dct(x);
% rng("default");
order = randperm(N, M);
z = z_full(order);
% get clipped threshold
epsilon = get_threshold_by_cr(CR,z);
signal = clip(z, epsilon);
[y, alpha, sigma] = linear_model_clip(signal,SNR);   

% initialize the input to the LMMSE part
zPri = zeros(N, 1); 
vZpri = 1;
% output variance of estimation of NLE, i.e., prediction of MSE
vzpris = zeros(1, maxIter);
% output variance of estimation of LE
vxpris = zeros(1, maxIter);
% the true mse each iteration
mses = zeros(1, maxIter);
for i=1:maxIter
   fprintf("iteration number %d \n", i);
   vzpris(i) = max(1e-7, vZpri);
%    [zPost,vZpost] = Z_APP_Clip(epsilon, y, zPri(order), vZpri, sigma );
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
%    xOrth(xOrth < 1e-6) = 0;
   mses(i) = norm(xOrth - x)^2 / length(x);
   vZpri = vXorth;
   if vZpri < 1e-7
       fprintf("iteration end at vXpri: %d, vZpri: %d, mse: %d \n", vXpri, vZpri, mses(i));
       mses(i+1:end) =mses(i);
       vzpris(i+1:end) = vZpri;
       vxpris(i+1:end) = vXpri;
       break;
   end
   zPri = dct(xOrth);
   fprintf("vXpri: %d, vZpri: %d, mse: %d \n", vXpri, vZpri, mses(i));
%    vZpri = mses(i);
end

xOrth(xOrth < 1e-3) = 0;
errorSections = 0;
recoverPositions = zeros(L,1);
for i = 0:L-1
    oneBlock = xOrth(i * B+1: (i+1) * B);
    maxValueOneBlockIndex = find(oneBlock == max(oneBlock));
    recoverPosition = i * B + maxValueOneBlockIndex;
    if(isempty(recoverPosition))
        errorSections = errorSections + 1;
        continue;
    end
    recoverPositions(i+1) = recoverPosition;
    if (recoverPositions(i+1) - positions(i+1)) > 1e-6
        errorSections = errorSections + 1;
    end
end
errorSectionRate = errorSections / L * 100;
fprintf("errorSectionRate: %%%.2f",  errorSectionRate);



figure
semilogy(1:maxIter, vzpris, "-x");
hold on
semilogy(1:maxIter, mses, '-^');
legend("Predicted MSE", "True MSE");
xlabel("iteration number");
ylabel("MSE");

figure
semilogx(vzpris, vxpris, '-bd');
hold on
semilogx([vzpris(2:end), vZpri], vxpris, '-x');
legend("de-clip", "de-modulation");
% hold on
% load se_mat.mat
% semilogx(10.^ interpolate_range, 10 .^interpolated_vz_orths_nle, '-');
% hold on
% semilogx(10.^ interpolate_range, 10 .^ interpolated_vz_orths_le, '--');
% xlim([5e-5, 1]);
% legend("de-clip", "de-modulation", "simulated de-modulation", "simulated de-clip");
xlabel("vz");
ylabel("vx");
end

