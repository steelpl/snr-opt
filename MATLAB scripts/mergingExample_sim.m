clear; clc
%% This is an example of merging 2 synthetic zero mean datasets (no time series)
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
p = 2; % number of datasets
eta = ones(p,1);
rng(123); % random seed

%%%%%%%%% scaling factor (a)%%%%%%%%%%%%%%
a = ones(p,1); % 1-vector 
% a = [0.5;1.5]; % non-1-vector (for test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Step 2: Synthetic data generation
% Signal and error: y and e
EeeT = eye(p);
Ey2 = 1;
ExxT = EeeT+Ey2*(a*a');
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
Ey2_est = Ey2*0.5; % estimated signal power
[EeeT_est,theta_est,rho2_est] = ECVest(ExxT); % modified SNRest for TC-like estimation
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