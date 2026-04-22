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
% Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021).
% Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
% IEEE Trans Geosci Remote Sens
%
% If you use the methods presented in the paper and/or this example,
% please cite this paper where appropriate.
% ============================================================

%% ---------- Step 1: Data parameters ----------
p = 5;                 % number of datasets
eta = ones(p,1);
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
MSE = @(u) diag(u'*ExxT*u - 2*Ey2*u'*a + Ey2);
R2  = @(u) diag(Ey2*((u'*(a*a')*u) ./ max(u'*ExxT*u, 1e-12)));

% Performance of individual datasets
MSE_ori = MSE(eye(p))';
R2_ori  = R2(eye(p))';

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

% Relation between WA and SNRopt (eq. 9 in the reference)
EeeT_inv_eta = EeeT\eta;
EeeT_inv_a   = EeeT\a;
s = diag(((eta'*EeeT_inv_eta)/(1/Ey2 + a'*EeeT_inv_a))*a); % us = s*uw

% Performance of merged data using true parameters
MSE_true = MSE([uw,us,ur,ue])';
R2_true  = R2([uw,us,ur,ue])';

% Display merged performance
disp([newline,'+ Metrics for merged data by ''true'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ', num2str(round(MSE_true,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ', num2str(round(R2_true,3))])

%% ---------- Step 4: Merging using estimated parameters ----------
% This section evaluates the practical case in which the required
% statistics are not known and must be estimated.

Ey2_est = Ey2*0.5;        % rough estimate of signal power
covx = ExxT;              % cov(x) = E[xxT] because E[x]=0

% Estimated quantities
[EeeT_est, theta_est, rho2_est] = NC(covx);
[N_est, a_est] = SNRest(ExxT, Ey2_est);

% Merging weights based on estimated parameters
uw_est = WA(EeeT_est);        % weighted average
us_est = SNRopt(N_est,a_est); % SNR-opt
ur_est = maxR(a_est,ExxT);    % maxR weight

% Relation between WA and SNRopt under estimated parameters
EeeT_inv_a_est = EeeT\a_est;
s_est = diag(((eta'*EeeT_inv_eta)/(1/Ey2_est + a_est'*EeeT_inv_a_est))*a_est); % us = s*uw

% Performance of merged data using estimated parameters
MSE_est = MSE([uw_est,us_est,ur_est,ue])';
R2_est  = R2([uw_est,us_est,ur_est,ue])';

% Display merged performance
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
legend('Original','Merged by true prm','Merged by est prm')
title('MSE')
set(gca,'xticklabel',method)
axis square

subplot(1,2,2)
bar(R2_results')
ylabel('R^2')
title('R^2')
set(gca,'xticklabel',method)
axis square