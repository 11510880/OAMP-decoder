function [errorSectionRate] = omap_clip_irr(B,L,CRs,lambdas,N,R,SNR,maxIter)
[x, positions] = generate_sparse_verctor(B,L);
M = L * log2(B) / R;
z_full = dct(x);
% rng("default");
order = randperm(N, M);
%% 
z = z_full(order);
% get clipped threshold
counter = 0;
y = [];
lambdas(lambdas<1e-4) = 0;
zero_indexes = find(lambdas==0);
lambdas(zero_indexes) = [];
CRs(zero_indexes) = [];
clip_indexes = zeros(length(CRs)+1,1);
epsilons = zeros(length(CRs), 1);
sigmas = zeros(length(CRs), 1);
for i=1:length(lambdas)
    CR = CRs(i);
    lambda = lambdas(i);
    n = max(floor(M*lambda), 0);
    z_one_cr = z(counter+1:counter+n);
    counter = counter+n;
    clip_indexes(i+1) = counter;
    if i==length(lambdas)
        z_one_cr = [z_one_cr;z(counter+1:end)];
        clip_indexes(end) = M;
    end
    epsilons(i) = get_threshold_by_cr(CR,z_one_cr);
    signal_one_cr = clip(z_one_cr, epsilons(i));
%     signal = [signal;signal_one_clip];
    [y_one_cr, alpha, sigmas(i)] = linear_model_clip(signal_one_cr,SNR);
    y = [y;y_one_cr];
end


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
   vzpris(i) = max(1e-7, vZpri);
   zPosts = [];
   vZposts = [];
   for clip_i=1:length(CRs)
      clip_range = clip_indexes(clip_i)+1:clip_indexes(clip_i+1);
      zPriResample = zPri(order);
      [zPost,vZpost] = declip(y(clip_range), zPriResample(clip_range), vZpri, sigmas(clip_i), epsilons(clip_i));
%       [zPost,vZpost] = Z_APP_Clip(epsilons(clip_i), y(clip_range), zPriResample(clip_range), vZpri, sigmas(clip_i));
      zPosts = [zPosts; zPost];
      vZposts = [vZposts; vZpost];
   end
   zPost = zPosts;
   vZpost = vZposts' * lambdas;
   % [zPost,vZpost] = declip(y, zPri(order), vZpri, sigma, epsilon);
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
    recoverPositions(i+1) = recoverPosition;
    if (recoverPositions(i+1) - positions(i+1)) > 1e-6
        errorSections = errorSections + 1;
    end
end
errorSectionRate = errorSections / L * 100;
fprintf("errorSectionRate: %%%.2f",  errorSectionRate);
end

