clear; clc

%% 1. Read Excel file and preprocess
all = readmatrix('irrigatgion_sample.xlsx', 'Range', 'B:E');
all = all(2:end,:);                 % remove header row

col_means    = mean(all);           % mean of each column
all_centered = all - col_means;     % column-wise centering

x = all_centered(:,1:end-1);        % predictors (centered)
p = size(x,2);                      % number of predictors
y = all_centered(:,end);            % response (centered)

xb= x+col_means(1:end-1);
yb= y+col_means(end);

%% 2. Calculate true parameters (a, N)
Ey2  = var(y,1);                    % signal power
n    = size(x,1);                   % Data length
cxy  = (x.'*y)/n;                   % Covariance b/w each x and y
a = cxy / Ey2;                      % scaling factor
e   = y - x.* a.';                  % error (per column)
EeeT = cov(e);                      % error covariance (E[ee^T])
N    = EeeT ./ Ey2;                 % error-to-signal ratio

%% 3. Estimate parameters through SNR-est (a_est, N_est)
Ey2_est = Ey2;                      % normally from long-term data, here assumed known                                    
ExxT = cov(x,1);                   % covariance of x (E[xx^T])
[N_est,a_est] = SNRest(ExxT, Ey2_est);

%% 4. Data merging
% 4.1 using true parameters
u      = SNRopt(N,a);               % optimal weights (true parameters)
ym     = x * u;                     % merged data with true params
bm     = col_means(1:end-1) * u;    % bias with u
ymb    = ym + bm;                   % adding bias to ym

% 4.2 using estimated parameters
u_est  = SNRopt(N_est,a_est);       % optimal weights (estimated parameters)
ym_est = x * u_est;                 % merged data with estimated params
bm_est = col_means(1:end-1) * u_est;% bias with u_est
ymb_est= ym_est + bm_est;           % adding bias ym_est

%% 5. Performance metrics (RMSE, R)
% 5.1 Using weights (identical results with 5.2)
% RMSE = @(u) sqrt(u'*ExxT*u - 2*Ey2*u'*a + Ey2);     % RMSE of estimator u'x vs y
% R = @(u) (sqrt(Ey2) * (u.'*a)) ./ sqrt(u.'*ExxT*u); % R b/w estimator u'x and y
% 
% RMSE_all = diag(RMSE([eye(p), u, u_est]))';    % RMSE for originals + merged
% R_all  = diag(R([eye(p), u, u_est]))';         % R  for originals + merged

% 5.2 Direct calculation
% 5.2.1 For debiased data
RMSE_all = sqrt(mean((y-[x,ym,ym_est]).^2)); 
R_all = corr(y,[x,ym,ym_est]); 

% 5.2.1 For biased data
RMSEb_all = mean((yb-[xb,ymb,ymb_est]).^2); 
Rb_all = corr(yb,[xb,ymb,ymb_est]); 

%% 6. Print results
varNames = [arrayfun(@(i) sprintf('Original_%d',i),1:p,'UniformOutput',false), ...
            {'SNRopt_true','SNRopt_est'}];

% 6.1 For debiased data (SNRopt results)
RMSE_fmt = round(RMSE_all(:),3);
R_fmt  = round(R_all(:),3);

T = table(RMSE_fmt, R_fmt, ...
          'VariableNames',{'RMSE','R'}, ...
          'RowNames',varNames);

fprintf('\n--- Merging results for debiased data (SNRopt results) ---\n');
disp(T)

% 6.2 For biased data (after adding bias => likely to be suboptimal)
RMSEb_fmt = round(RMSEb_all(:),3);
Rb_fmt  = round(Rb_all(:),3);

Tb = table(RMSEb_fmt, Rb_fmt, ...
          'VariableNames',{'RMSE','R'}, ...
          'RowNames',varNames);

fprintf('\n--- Merging results for biased data (after adding bias => likely to be suboptimal) ---\n');
disp(Tb)