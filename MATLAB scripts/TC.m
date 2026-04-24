function [EeeT_est, theta_est, rho2_est] = TC(covx)
% ============================================================
% Triple Collocation (TC)
% Closed-form estimation assuming mutually uncorrelated errors
%
% INPUT
%   covx : 3x3 covariance matrix of three datasets
%
% OUTPUT
%   EeeT_est  : 3x3 diagonal matrix of estimated error variances
%   theta_est : 3x1 relative scaling vector (theta(1) = 1)
%   rho2_est  : 3x1 squared correlation between each dataset and the latent signal
%
% METHOD SUMMARY
%   Model: Q = theta * theta.' * Var(y) + EeeT, with EeeT diagonal.
%   Using the three off-diagonal covariances Q12,Q13,Q23, closed-form
%   expressions yield estimates of error variances and the common
%   signal variance. Scaling is fixed by setting theta(1)=1.
%
% VALIDITY AND NOTES
%   - Function only supports p = 3.
%   - Requires nonzero off-diagonal covariances; otherwise returns NaNs.
%   - If any estimated error variance is negative or non-finite, the
%     entire solution is invalid and NaNs are returned.
%   - rho2 = 1 - Eii / Var(y) is clipped for tiny numerical violations.
% ============================================================

%% ---------- Input check ----------
if nargin < 1
    error('TC requires one input argument: covx');
end
if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end

p = size(covx,1);
if p ~= 3
    error('TC is only valid for p = 3.');
end

%% ---------- Setup ----------
% Symmetrize to remove tiny asymmetry and extract entries
Q = (covx + covx.') / 2;

Q11 = Q(1,1); Q22 = Q(2,2); Q33 = Q(3,3);
Q12 = Q(1,2); Q13 = Q(1,3); Q23 = Q(2,3);

tol = 1e-12;

%% ---------- Validity check for denominators ----------
% TC requires nonzero cross-covariances
if abs(Q12) < tol || abs(Q13) < tol || abs(Q23) < tol
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Error variance estimation ----------
% Closed-form TC estimates assuming zero error cross-covariances
E1 = Q11 - (Q12 * Q13) / Q23;
E2 = Q22 - (Q12 * Q23) / Q13;
E3 = Q33 - (Q13 * Q23) / Q12;
Evec = [E1; E2; E3];

%% ---------- Strict validity condition ----------
% Any negative or non-finite error variance invalidates the solution
if any(Evec < 0) || any(~isfinite(Evec))
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

% Diagonal error covariance matrix
EeeT_est = diag(Evec);

%% ---------- Signal variance estimation ----------
% Use cyclic relation; choose the (1,2)-(1,3)/(2,3) formulation
Ey2_est = (Q12 * Q13) / Q23;

if ~(isfinite(Ey2_est) && Ey2_est > 0)
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Scaling factor estimation ----------
% Normalize to dataset 1 (theta_1 = 1)
theta_est = [1;
             Q23 / Q13;
             Q23 / Q12];

%% ---------- Squared correlation estimation ----------
rho2_est = 1 - Evec / Ey2_est;

% Clip tiny numerical violations to [0,1]
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;