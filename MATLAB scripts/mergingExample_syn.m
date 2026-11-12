clear; clc
%% This is an example of merging 3 synthetic zero mean datasets (time series)
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
%
%% Step 1: Data parameters
p = 3; % number of datasets
n = 100; % data length
rng(123); % random seed

%%%%%%%%% scaling factor (a)%%%%%%%%%%%%%%
% a = ones(p,1); % 1-vector
a = rand(p,1); % non-1-vector (for test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecc = 0.3; % error cross-correlation [0,1]
SNRdB = 0.1; % SNR in dB

%% Step 2: Synthetic data generation
% Signal and error: y and e
[y,e] = dataGEN(n,p,ecc,SNRdB);
% SNRdB_actual = snr(y, e(:,1)); % for checking produced SNR

% Observation: x = a*signal + error
x = y*a' + e;

% Signal power and covariance matrices of e and x
Ey2 = var(y); % signal power, E[y^2] = var(y) because E[y]=0
EeeT = cov(e); % E[eeT] = cov(e) because E[e]=0
ExxT = cov(x); % E[xxT] = cov(x) because E[x]=0
N = EeeT./Ey2; % error-to-signal ratio

% Functions of MSE and R2 by linear combination using u
MSE = @(u) sum(u.*(ExxT*u),1) - 2*Ey2*(a'*u) + Ey2;
R2 = @(u) Ey2*(a'*u).^2 ./ sum(u.*(ExxT*u),1);

% MSE of observations
MSE_ori = MSE(eye(p));

% Pearson correlation of observations
R2_ori = R2(eye(p));

% Printing metrics
disp('+ Metrics for original data')
disp([' * MSE for x: ',num2str(round(MSE_ori,3))])
disp([' * R2 for x: ',num2str(round(R2_ori,3))])

%% Step 3: Merging using true parameters

% Merging weights and centered series
uw = WA(EeeT);  % weighted average
us = SNRopt(N,a); % SNR-opt
ur = maxR(a,ExxT); % maxR weight
ue = 1/p*ones(p,1); % equal weight
[yw,ys,yr,ye] = deal(x*uw,x*us,x*ur,x*ue);

% MSE of data merged by true parameters
MSE_true = MSE([uw,us,ur,ue]);

% Pearson correlation of data merged by true parameters
R2_true = R2([uw,us,ur,ue]);

% Printing metrics
disp([newline,'+ Metrics for merged data by ''true'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ',num2str(round(MSE_true,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ',num2str(round(R2_true,3))])

%% Step 4: Merging using estimated parameters
% Estimation of merging statistics
Ey2_est = 0.5 * Ey2; % deliberately rough external-scale proxy in this simulation
[N_est, a_est, snr_est, fit] = SNRest(ExxT, Ey2_est);

% Reuse the same fitted error and observation covariances for all weights.
uw_est = WA(fit.errorCov,fit.info.covarianceScale);
us_est = SNRopt(N_est,a_est);
ur_est = maxR(fit.theta,fit.covariance);
[yw_est,ys_est,yr_est] = deal(x*uw_est,x*us_est,x*ur_est);

MSE_est = MSE([uw_est,us_est,ur_est,ue]);
R2_est = R2([uw_est,us_est,ur_est,ue]);

fprintf('\nNC min eig(E): %.6g\n', fit.info.minErrorCovEigenvalue);
fprintf('NC linear SNR by product: %s\n', mat2str(snr_est',4));

disp([newline,'+ Metrics for merged data by ''estimated'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ',num2str(round(MSE_est,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ',num2str(round(R2_est,3))])

%% Step 5: Plotting merging results
MSE_results = [mean(MSE_ori)*ones(1,4);MSE_true;MSE_est];
R2_results = [mean(R2_ori)*ones(1,4);R2_true;R2_est];
method = {'WA','SNRopt','maxR', 'EW'}; % weighted average, SNRopt, maxR, Equal Weight

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
