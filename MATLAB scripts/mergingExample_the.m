clear; clc
%% ============================================================
% Example: merging synthetic zero-mean datasets (no time series)
%
% NOTE
% (1) If a = 1, WA and SNRopt can be compared more fairly because WA does
%     not explicitly account for scaling factors.
% (2) If a is non-uniform, SNRopt and maxR can exploit scaling differences,
%     whereas WA cannot.
% (3) Ey2 should be reasonably estimated for SNRopt in practical settings
%     (e.g. reanalysis or climatology). Here, 0.5 of the true Ey2 is used
%     as an intentionally rough estimate.
%
% REFERENCE
% For more details, see:
%
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Transactions on Geoscience and Remote Sensing, 60, 1–15.
%   https://doi.org/10.1109/TGRS.2021.3107028
%
%   Pham, H. T., & Kim, S. (2027). Journal of Hydrology, 682, 136749.
%   https://doi.org/10.1016/j.jhydrol.2026.136749
%
% If you use the methods presented in the paper and/or this example,
% please cite this paper where appropriate.
% ============================================================

%% ---------- Step 1: Data parameters ----------
p = 5;                 % number of datasets
rng(123);              % random seed

%%%%%%%%% scaling factor (a)%%%%%%%%%%%%%%
% a = ones(p,1);       % 1-vector
a = rand(p,1);         % non-1-vector (for test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecc = 0.3;             % error cross-correlation [0,1]
SNRdB = 0.1;           % SNR in dB

%% ---------- Step 2: Synthetic data generation ----------
% Construct the synthetic covariance structure:
%   ExxT = EeeT + Ey2 * (a a')
% where EeeT is the error covariance matrix and Ey2 is the signal power.

EeeT = EeeTGEN(p,ecc);             % error covariance matrix
SNR = 10^(SNRdB/10);               % SNR conversion
Ey2 = mean(diag(EeeT)*SNR);        % signal power based on SNR and EeeT
ExxT = EeeT + Ey2*(a*a');          % E[xxT] = cov(x) because E[x]=0
N = EeeT./Ey2;                     % error-to-signal ratio matrix

% Functions of MSE and R2 by linear combination using u.
% These definitions support both a single weight vector (p x 1) and a
% matrix of weight vectors (p x m), returning one value per column.
MSE = @(u) sum(u.*(ExxT*u),1) - 2*Ey2*(a'*u) + Ey2;
R2 = @(u) Ey2*(a'*u).^2 ./ sum(u.*(ExxT*u),1);

% Performance of individual datasets
MSE_ori = MSE(eye(p));
R2_ori  = R2(eye(p));

% Display original performance
disp('+ Metrics for original data')
disp([' * MSE for x: ', num2str(round(MSE_ori,3))])
disp([' * R2 for x: ', num2str(round(R2_ori,3))])

%% ---------- Step 3: Merging using true parameters ----------
% This section evaluates the ideal performance when the true statistics
% (EeeT, Ey2, N, and a) are fully known.

uw = WA(EeeT);            % weighted average
us = SNRopt(N,a);         % SNR-opt
ur = maxR(a,ExxT);        % maxR weight
ue = ones(p,1)/p;         % equal weight

% Performance of merged data using true parameters
MSE_true = MSE([uw,us,ur,ue]);
R2_true  = R2([uw,us,ur,ue]);

% Display merged performance
disp([newline,'+ Metrics for merged data by ''true'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ', num2str(round(MSE_true,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ', num2str(round(R2_true,3))])

%% ---------- Step 4: Merging using estimated parameters ----------
% This section evaluates the practical case in which the required
% statistics are not known and must be estimated.

Ey2_est = 0.5 * Ey2; % deliberately rough external-scale proxy in this simulation
[N_est, a_est, snr_est, fit] = SNRest(ExxT, Ey2_est);

% Reuse the same fitted error and observation covariances for all weights.
uw_est = WA(fit.errorCov,fit.info.covarianceScale);
us_est = SNRopt(N_est,a_est);
ur_est = maxR(fit.theta,fit.covariance);

MSE_est = MSE([uw_est,us_est,ur_est,ue]);
R2_est = R2([uw_est,us_est,ur_est,ue]);

fprintf('\nNC min eig(E): %.6g\n', fit.info.minErrorCovEigenvalue);
fprintf('NC linear SNR by product: %s\n', mat2str(snr_est',4));

disp([newline,'+ Metrics for merged data by ''estimated'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ', num2str(round(MSE_est,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ', num2str(round(R2_est,3))])

%% ---------- Step 5: Plotting merging results ----------
% Compare the average performance of original datasets against the merged
% products using true and estimated parameters.

MSE_results = [mean(MSE_ori)*ones(1,4); MSE_true; MSE_est];
R2_results  = [mean(R2_ori)*ones(1,4); R2_true; R2_est];

method = {'WA','SNRopt','maxR','EW'}; % Weighted Average, SNRopt, maxR, Equal Weight

figure('Units','normalized','Position',[0 0 1 1])

subplot(1,2,1)
bar(MSE_results')
ylabel('MSE')
legend('Original','Merged by true prm','Merged using NC')
title('MSE')
set(gca,'xticklabel',method)
axis square

subplot(1,2,2)
bar(R2_results')
ylabel('R^2')
title('R^2')
set(gca,'xticklabel',method)
axis square
