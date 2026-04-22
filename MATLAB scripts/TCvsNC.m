clear; clc
%% ECC sweep analysis: TC vs NC

p = 3;
rng(123);

ecc_list = [0:0.1:0.9,0.99];
nE = length(ecc_list);

SNRdB = 0.1;
SNR = 10^(SNRdB/10);

nMC = 1000;
tol = 1e-10;

%% Common random numbers across ecc
A_all    = rand(p,nMC);
varE_all = rand(p,nMC);

%% Storage
% (1) TC - NC difference (diag)
diff_med = nan(nE,1);
diff_q1  = nan(nE,1);
diff_q3  = nan(nE,1);

% (2) MSE vs TRUE
mse_tc_med = nan(nE,1);
mse_tc_q1  = nan(nE,1);
mse_tc_q3  = nan(nE,1);

mse_nc_med = nan(nE,1);
mse_nc_q1  = nan(nE,1);
mse_nc_q3  = nan(nE,1);

% (3) failure rate
tc_fail_rate = nan(nE,1);
nc_fail_rate = nan(nE,1);

%% Loop over ecc
for ie = 1:nE
    
    ecc = ecc_list(ie);
    
    % Preallocation
    diff_store  = nan(nMC,p);   % diagonal-wise difference
    mse_tc_all  = nan(nMC,1);
    mse_nc_all  = nan(nMC,1);
    
    tc_fail = 0;
    nc_fail = 0;
    
    for k = 1:nMC
        
        % --- generate system (same realization across ecc) ---
        a = A_all(:,k);
        var_e = varE_all(:,k);
        
        EeeT_true = EeeTGEN(p,ecc,var_e);
        
        Ey2 = mean(diag(EeeT_true)*SNR);
        ExxT = EeeT_true + Ey2*(a*a');
        
        % --- TC / NC ---
        [EeeT_TC, ~, ~] = TC(ExxT);
        [EeeT_NC, ~, ~] = NC(ExxT);
        
        % --- TC failure ---
        tc_invalid = any(isnan(EeeT_TC(:))) || any(diag(EeeT_TC) < 0);
        if tc_invalid
            tc_fail = tc_fail + 1;
        end
        
        % --- NC failure (tolerance-based) ---
        nc_invalid = any(isnan(EeeT_NC(:))) || any(diag(EeeT_NC) < -tol);
        if nc_invalid
            nc_fail = nc_fail + 1;
        end
        
        % --- diagonal ---
        E_true = diag(EeeT_true);
        
        % --- MSE vs TRUE (independent validity) ---
        if ~tc_invalid
            E_TC = diag(EeeT_TC);
            mse_tc_all(k) = mean((E_TC - E_true).^2);
        end
        
        if ~nc_invalid
            E_NC = diag(EeeT_NC);
            mse_nc_all(k) = mean((E_NC - E_true).^2);
        end
        
        % --- difference (TC - NC): only when both valid ---
        if ~tc_invalid && ~nc_invalid
            E_TC = diag(EeeT_TC);
            E_NC = diag(EeeT_NC);
            diff_store(k,:) = E_TC - E_NC;
        end
        
    end
    
    % ===== statistics =====
    
    % diff: flatten valid values only
    diff_all = diff_store(~isnan(diff_store));
    if ~isempty(diff_all)
        diff_med(ie) = median(diff_all);
        diff_q1(ie)  = prctile(diff_all,25);
        diff_q3(ie)  = prctile(diff_all,75);
    end
    
    % MSE TC
    mse_tc_valid = mse_tc_all(~isnan(mse_tc_all));
    if ~isempty(mse_tc_valid)
        mse_tc_med(ie) = median(mse_tc_valid);
        mse_tc_q1(ie)  = prctile(mse_tc_valid,25);
        mse_tc_q3(ie)  = prctile(mse_tc_valid,75);
    end
    
    % MSE NC
    mse_nc_valid = mse_nc_all(~isnan(mse_nc_all));
    if ~isempty(mse_nc_valid)
        mse_nc_med(ie) = median(mse_nc_valid);
        mse_nc_q1(ie)  = prctile(mse_nc_valid,25);
        mse_nc_q3(ie)  = prctile(mse_nc_valid,75);
    end
    
    % failure rate (%)
    tc_fail_rate(ie) = tc_fail / nMC * 100;
    nc_fail_rate(ie) = nc_fail / nMC * 100;
    
end

%% =========================
%% Visualization
%% =========================

figure

%% (1,3,1) Difference (TC - NC)
subplot(1,3,1)
hold on

fill([ecc_list fliplr(ecc_list)], ...
     [diff_q1'*1e3 fliplr(diff_q3'*1e3)], ...
     [0.7 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.5)

plot(ecc_list, diff_med*1e3, 'b-', 'LineWidth',2)

ylim([-10 10])
xlabel('ecc')
ylabel('TC - NC (×10^{-3})')
title('Difference (TC - NC)')
grid on
axis square 

%% (1,3,2) MSE vs TRUE
subplot(1,3,2)
hold on

fill([ecc_list fliplr(ecc_list)], ...
     [mse_tc_q1' fliplr(mse_tc_q3')], ...
     [1 0.7 0.7], 'EdgeColor','none', 'FaceAlpha',0.5)
plot(ecc_list, mse_tc_med, 'r-', 'LineWidth',2)

fill([ecc_list fliplr(ecc_list)], ...
     [mse_nc_q1' fliplr(mse_nc_q3')], ...
     [0.7 0.8 1], 'EdgeColor','none', 'FaceAlpha',0.5)
plot(ecc_list, mse_nc_med, 'b-', 'LineWidth',2)

xlabel('ecc')
ylabel('Error Variance MSE')
title('MSE of Estimated Error Variance (TC vs NC)')
legend('TC IQR','TC median','NC IQR','NC median','Location','northwest')
grid on
axis square 

%% (1,3,3) Failure rate
subplot(1,3,3)
hold on

plot(ecc_list, tc_fail_rate, 'r-o', 'LineWidth',2)
plot(ecc_list, nc_fail_rate, 'b-o', 'LineWidth',2)

xlabel('ecc')
ylabel('Failure rate (%)')
title('Failure Rate')
legend('TC','NC','Location','northwest')
grid on
axis square 