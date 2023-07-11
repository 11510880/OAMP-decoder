clear all;
clc;
% constants
B = 64;
L = 2048;
CR = -3;
SNRs = 7:0.25:9;
SNRs_irr = 5:0.25:7;
R = 1;
maxIter = 50;
N = B * L;
rng("shuffle");
simulationTimes = 1;
SER_rclip = zeros(length(SNRs), simulationTimes);
SER_irclip = zeros(length(SNRs_irr), simulationTimes);
numLevels = 64;
for simuTime=1:simulationTimes
    SERs = zeros(length(SNRs), 1);
    for i=1:length(SNRs)
        SNR = SNRs(i);
%         errorSectionRate = oamp_clip(B,L,CR,N,R,SNR,maxIter);
        errorSectionRate=oamp_quant(B,L,numLevels,N,R,SNR,maxIter);
        SERs(i) = errorSectionRate;
    end
    
    SERs(SERs<1e-5) = 1e-5;
    SER_rclip(:,simuTime) = SERs;
    
    %% irregular clip
    load("clippedRatios.mat");
    load("lambdas.mat")
    
    for i=1:length(SNRs_irr)
        SNR = SNRs_irr(i);
        errorSectionRate = omap_clip_irr(B,L,clippedRatios,lambdas,N,R,SNR,maxIter);
        SERs(i) = errorSectionRate;
    end
    
    SERs(SERs<1e-5) = 1e-5;
    SER_irclip(:,simuTime) = SERs;
end
figure
semilogy(SNRs, mean(SER_rclip,2), '-x');
hold on
semilogy(SNRs_irr, mean(SER_irclip,2), '-^');
xlabel("SNR(dB)");
ylabel("Section Error Rate(SER)");
legend("Regular clip", "Irregular clip");



