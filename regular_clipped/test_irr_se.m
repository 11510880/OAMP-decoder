% test state evolution of irregular clipping
%% constants
rng('default');
SNR = 6.4      ;
B = 64;
L = 2048;
R = 1.0;
% clippedRatios = [-300,-12,-10,-8];
clippedRatios = [-300,-30, -22, -20, -19, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 300];


[lambdas, le_variance_x_optimal, nle_variance_x] = state_evolutions(SNR, clippedRatios, B, L, R);
interpolate_range = -5:0.01:0;
variance_z = 10.^ interpolate_range;
semilogx(variance_z(1:486), le_variance_x_optimal); hold on; semilogx(variance_z(1:486), nle_variance_x); legend("le", "nle");