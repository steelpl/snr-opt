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
%% Step 1: Data parameters
p = 3; % number of datasets
eta = ones(p,1);
rng(123); % random seed

%%%%%%%%% scaling factor (a)%%%%%%%%%%%%%%
% a = ones(p,1); % 1-vector 
a = rand(p,1); % non-1-vector (for test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecc = 0.3; % error cross-correlation [0,1]
SNRdB = 0.1; % SNR in dB

%% Step 2: Synthetic data generation
% error covariance matrix, signal power, N
EeeT = EeeTGEN(p,ecc); % error covariance matrix
SNR = 10^(SNRdB/10); % SNR conversion
Ey2 = mean(diag(EeeT)*SNR); % Signal power based on given SNR and EeeT
ExxT = EeeT+Ey2*(a*a'); % E[xxT] = cov(x) because E[x]=0
N = EeeT./Ey2; % error-to-signal ratio

% Functions of MSE and R2 by linear combination using u
MSE = @(u) u'*ExxT*u - 2*Ey2*u'*a + Ey2;
R2 = @(u) Ey2*((u'*(a*a')*u)./(u'*ExxT*u));

% MSE of observations
MSE_ori = diag(MSE(eye(p)))';

% Squared Pearson correaltion of observations
R2_ori = diag(R2(eye(p)))';

% Printing metrics
disp('+ Metrics for original data')
disp([' * MSE for x: ',num2str(round(MSE_ori,3))])
disp([' * R2 for x: ',num2str(round(R2_ori,3))])

%% Step 3: Merging using true parameters

% weight of three methods
uw = WA(EeeT); % weighted average
us = SNRopt(N,a); % SNR-opt
ur = maxR(a,ExxT);% maxR weight
ue = 1/p*ones(p,1);% equal weight 
s = diag(((eta'*(EeeT\eta))/(1/Ey2+a'*(EeeT\a)))*a); % us = s*uw (eq. 9)

% MSE of merged data
MSE_true = diag(MSE([uw,us,ur,ue]))';

% Squared Pearson correaltion of merged data
R2_true = diag(R2([uw,us,ur,ue]))';

% Printing metrics
disp([newline,'+ Metrics for merged data by ''true'' parameters'])
disp([' * MSE for WA, SNRopt, maxR, EW: ',num2str(round(MSE_true,3))])
disp([' * R2 for WA, SNRopt, maxR, EW: ',num2str(round(R2_true,3))])

%% Step 4: Merging using estimated parameters
% merging statistics
Ey2_est = Ey2*0.5; % roughly estimated signal power (e.g., reanalysis, var(mean(x,2)). Here, 0.5 of the true Ey2 is arbitrarily selected.
covx=ExxT; % cov(x) = E[xxT] because E[x]=0
[EeeT_est,theta_est,rho2_est] = ECVest(covx); % modified SNRest for TC-like estimation
[N_est,a_est] = SNRest(ExxT, Ey2_est); % SNR-est

% weight of three methods
uw_est = WA(EeeT_est); % weighted average
us_est = SNRopt(N_est,a_est); % SNR-opt
ur_est = maxR(a_est,ExxT);% maxR weight
s_est = diag(((eta'*(EeeT\eta))/(1/Ey2+a_est'*(EeeT\a_est)))*a_est); % us = s*uw (eq. 9)

% RMSE of merged data
MSE_est = diag(MSE([uw_est,us_est,ur_est,ue]))';

% Pearson correaltion of merged data
R2_est = diag(R2([uw_est,us_est,ur_est,ue]))';

% Printing metrics
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