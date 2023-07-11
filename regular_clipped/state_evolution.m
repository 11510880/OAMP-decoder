clear all
clc
%% constants
rng('default');
SNR = 2.0;
B = 64;
L = 2048;
R = 0.5;
N = B * L;
M = L * log2(B) / R;
CR = -13;
interpolate_range = -5:0.01:0;
simulation_times = 5;
%% LE simulation
% vz from 0 to 60dB
rho = 0:60;
% 10log10(1/vz) = rho
vz_pris = 1./(10.^(rho/10));
vz_orths_history = zeros(length(vz_pris), simulation_times);
for simu_time=1:simulation_times
    fprintf("repetion %d \n", simu_time);
    vz_posts = zeros(length(vz_pris), 1);
    vz_orths = zeros(length(vz_pris), 1);
    
    for i = 1:length(vz_pris)
        vz_pri = vz_pris(i);
        z_pri = sqrt(1-vz_pri) * randn(N, 1);
        n_pri = sqrt(vz_pri) * randn(N, 1);
        % here p(z) ~ N(zpri, npri)
        z = z_pri + n_pri;
        order = randperm(N, M);
        z_resample = z(order);
        % get clipped threshold
        epsilon = get_threshold_by_cr(CR, z_resample);
%         epsilon = 0.0496;
        signal = clip(z_resample, epsilon);
        [y, alpha, sigma] = linear_model_clip(signal,SNR);
%         [z_post,vz_post] = declip(y, z_pri(order), vz_pri, sigma, epsilon);
        [z_post,vz_post] = Z_APP_Clip( epsilon, y, z_pri(order), vz_pri, sigma );
        [vz_orth, z_orth] = orthogonalization(vz_post, vz_pri, z_post, z_pri, order);
        vz_posts(i) = vz_post;
        vz_orths(i) = vz_orth;
    end
    vz_orths_history(:,simu_time) = vz_orths;
end
vz_orths = mean(vz_orths_history, 2);
vz_pris_log = log10(vz_pris);
vz_orths_log = log10(vz_orths);
% interpolated_vz_orths_le = interp1(vz_pris_log, log10(vz_orths_history(:,1)), interpolate_range);
interpolated_vz_orths_le = interp1(vz_pris_log, vz_orths_log, interpolate_range);

%% NLE simulation
vx_pris_log = [[-1.2:0.01:0-0.01]'; [0:0.05:1-0.05]'; [1:0.1:6]'];
vx_pris_log = vx_pris_log(end:-1:1);
vx_pris = 10 .^ vx_pris_log;
vx_orths_history = zeros(length(vx_pris),  simulation_times);
mse_history = zeros(length(vx_pris),  simulation_times);
for simu_time=1:simulation_times
    fprintf("repetion %d \n", simu_time);
    x = generate_sparse_verctor(B,L);
    vx_posts = zeros(length(vx_pris), 1);
    vx_orths = zeros(length(vx_pris), 1);
    mses = zeros(length(vx_pris), 1);
    for i=1:length(vx_pris)
        vx_pri = vx_pris(i);
        x_pri = x + sqrt(vx_pri) * randn(N, 1);
        [x_post,vx_posts(i)] = sr_demodulation(B, L, vx_pri, x_pri);
        % xPost N * 1,  xPri, N * 1
        [vx_orths(i), x_orth] = orthogonalization(vx_posts(i),vx_pri, x_post, x_pri);
        mses(i) = mean((x_orth - x).^2);
    end
    vx_orths_history(:,simu_time) = vx_orths;
    mse_history(:, simu_time) = mses;
end
vx_orths = mean(vx_orths_history, 2);
mse_average = mean(mse_history, 2);
zero_index = find(vx_orths == 0);
vx_orths(zero_index) = [];
vx_pris_log(zero_index) = [];
vx_orths_log = log10(vx_orths);
interpolated_vz_orths_nle = interp1(vx_orths_log(50:end), vx_pris_log(50:end), interpolate_range);
figure
plot(vx_orths, '-x');
hold on;
plot(mse_average, '-^');
legend("predicted MSE", "true MSE");

%%
figure
semilogx(10.^ interpolate_range, 10 .^interpolated_vz_orths_nle, 'LineWidth', 1.5);
hold on
semilogx(10.^ interpolate_range, 10 .^ interpolated_vz_orths_le, 'LineWidth', 1.5);
xlabel("vz");
ylabel("vx");
legend("NLE", "LE");
save("se_mat.mat", "interpolate_range", "interpolated_vz_orths_nle", "interpolated_vz_orths_le");







