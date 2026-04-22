clear; clc
%% This is an example of merging 3 synthetic zero mean datasets (no time series)
% NOTE:
% (1) This example has a scaling factor (a) equal to 1, which is for fair 
% comparison of the two merging methods (WA and SNRopt)
% (2) Since WA does not consider 'a' in its merging process, a fair
% comparison with SNRopt is not available unless the scaling factors are 
% the same sa each other. On contrary, SNRopt allows arbitrary 'a' for 
% weight computation. Any non-1 scaling factor can be tested in this 
% example to see how the merging results are different.
% (3) Ey2 (signal power) shoud be reasonably estimated for SNRopt 
% (e.g. reanalysis or climatology). Here, 0.5 of true Ey2 is tested.
%
% REFERENCE
% For more details, see:
%
% Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021). 
% Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
% IEEE Trans Geosci Remote Sens
%
% If you use the methods presented in the paper and/or this example, 
% please cite this paper where appropriate.
%
clear; clc

%% Step 1: Data parameters
p = 3;
eta = ones(p,1);
rng(123);

% scaling factor
% a = ones(p,1);
a = rand(p,1);

ecc = 0.3;
SNRdB = 0.1;

%% Step 2: Synthetic data generation
EeeT = EeeTGEN(p,ecc);
SNR = 10^(SNRdB/10);
Ey2 = mean(diag(EeeT)*SNR);

ExxT = EeeT + Ey2*(a*a');
N = EeeT./Ey2;

% === More robust function definitions ===
MSE = @(u) diag(u'*ExxT*u - 2*Ey2*u'*a + Ey2);
R2  = @(u) diag(Ey2*((u'*(a*a')*u)./(u'*ExxT*u)));

% === Original data metrics (clean + stable) ===
U_eye = eye(p);

MSE_ori = MSE(U_eye)';
R2_ori  = R2(U_eye)';

disp('+ Metrics for original data')
disp([' * MSE for x: ', num2str(round(MSE_ori,3))])
disp([' * R2 for x: ', num2str(round(R2_ori,3))])

%% Step 3: Merging using true parameters

uw = WA(EeeT);
us = SNRopt(N,a);
ur = maxR(a,ExxT);
ue = ones(p,1)/p;

% (kept but numerically stabilized)
den = (1/Ey2 + a'*(EeeT\a));
num = (eta'*(EeeT\eta));
s = diag((num/den)*a);

% === vectorized evaluation ===
U_true = [uw, us, ur, ue];

MSE_true = MSE(U_true)';
R2_true  = R2(U_true)';

disp([newline,'+ Metrics for merged data by ''true'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ', num2str(round(MSE_true,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ', num2str(round(R2_true,3))])

%% Step 4: Merging using estimated parameters

Ey2_est = Ey2 * 0.5;
covx = ExxT;

[EeeT_est, theta_est, rho2_est] = NC(covx);
% [EeeT_est,theta_est,rho2_est] = ECVest(covx);

[N_est,a_est] = SNRest(ExxT, Ey2_est);

uw_est = WA(EeeT_est);
us_est = SNRopt(N_est,a_est);
ur_est = maxR(a_est,ExxT);

% (same stabilization)
den_est = (1/Ey2 + a_est'*(EeeT\a_est));
s_est = diag((num/den_est)*a_est);

U_est = [uw_est, us_est, ur_est, ue];

MSE_est = MSE(U_est)';
R2_est  = R2(U_est)';

disp([newline,'+ Metrics for merged data by ''estimated'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ', num2str(round(MSE_est,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ', num2str(round(R2_est,3))])

%% Step 5: Plotting merging results (unchanged logic, safer)

MSE_results = [mean(MSE_ori)*ones(1,4); MSE_true; MSE_est];
R2_results  = [mean(R2_ori)*ones(1,4); R2_true; R2_est];

method = {'WA','SNRopt','maxR','EW'};

figure('Units','normalized','Position',[0 0 1 1])

subplot(1,2,1)
bar(MSE_results')
ylabel('MSE')
legend('Original','Merged by true prm','Merged by est prm','Location','best')
title('MSE')
set(gca,'xticklabel',method)
axis square
grid on

subplot(1,2,2)
bar(R2_results')
ylabel('R^2')
title('R^2')
set(gca,'xticklabel',method)
axis square
grid on