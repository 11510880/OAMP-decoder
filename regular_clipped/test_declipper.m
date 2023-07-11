% % test declip and Z_APP_Clip
% load('declip.mat');
% [zPost,vPost, secondMoment, marginalLikelihood] = declip(y, zPri(order), vZpri, sigma, epsilon);
% [zPost1,vZpost1,integral_x2, p_y_gx] = Z_APP_Clip(epsilon, y, zPri(order), vZpri, sigma );
clear all;
clc;
% constants
B = 64;
L = 2048;
CR = -3;
SNR = 7;
R = 1;
maxIter = 50;
N = B * L;
rng("shuffle");
numLevels = 48;
[x, positions] = generate_sparse_verctor(B,L);
M = L * log2(B) / R;
z_full = dct(x);
vZpris = 0.1:0.1:2.0;
vZposts = zeros(length(vZpris),1);
vars = zeros(length(vZpris),1);
%% clip and declip
% get clipped threshold
% for i=1:length(vZpris)
%     vZpri = vZpris(i);
%     noise = sqrt(vZpri) * randn(length(z_full), 1);
%     z_gauss = z_full + noise;
%     epsilon = get_threshold_by_cr(CR,z_gauss);
%     signal = clip(z_gauss, epsilon);
%     [y, alpha, sigma] = linear_model_clip(signal,SNR);
%     [zPost,vZpost] = declip(y, z_full, vZpri, sigma, epsilon);
%     vZposts(i) = vZpost;
%     vars(i) = var(zPost - z_gauss);
% end
%% quant and dequant
for i=1:length(vZpris)
    vZpri = vZpris(i);
    noise = sqrt(vZpri) * randn(length(z_full), 1);
    z_gauss = z_full + noise;
    quant_max = max(z_gauss);
%     quant_max = prctile(z_gauss, 99.99);
    [z_full_quant, quant_levels,quant_bounds] = uniform_quantizer(z_gauss, numLevels, quant_max);
    [y, alpha, sigma] = linear_model_clip(z_full_quant, SNR);
    [zPost,vZpost] = de_quant(y, z_gauss, vZpri, sigma, quant_levels,quant_bounds);
    vZposts(i) = vZpost;
    vars(i) = var(zPost - z_gauss);
end

plot(vZpris, vars, '-*');
hold on;
plot(vZpris, vZposts,'-^');
legend("var(zPost - zTrue)", "vZposts");