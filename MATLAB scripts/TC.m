function [EeeT_est, theta_est, rho2_est] = TC(covx)
% ============================================================
% Triple Collocation (TC)
% Closed-form estimation under mutually uncorrelated errors
%
% INPUT
%   covx : covariance matrix (3 x 3)
%
% OUTPUT
%   EeeT_est  : estimated error covariance matrix (3 x 3, diagonal)
%   theta_est : estimated scaling vector (3 x 1)
%   rho2_est  : squared correlation with latent signal (3 x 1)
%
% METHOD
%   Classical TC assumes
%       Q = theta*theta' + EeeT
%   with zero error cross-covariances, i.e.,
%       offdiag(EeeT) = 0
%   and derives closed-form estimates from the covariance relations.
%
% NOTES
%   - TC is only valid for p = 3.
%   - If any estimated error variance becomes negative, the entire
%     solution is considered invalid and NaN is returned.
%   - Unlike NC, TC does not enforce feasibility through projection.
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
% Symmetrize Q to suppress tiny numerical asymmetry in the input.
Q = (covx + covx') / 2;

% Extract covariance entries for readability.
Q11 = Q(1,1); Q22 = Q(2,2); Q33 = Q(3,3);
Q12 = Q(1,2); Q13 = Q(1,3); Q23 = Q(2,3);

tol = 1e-12;

%% ---------- Validity check for denominators ----------
% TC formulas rely on nonzero cross-covariances.
% If any denominator is numerically zero, the solution is invalid.
if abs(Q12) < tol || abs(Q13) < tol || abs(Q23) < tol
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Error variance estimation ----------
% Classical TC error variances under the assumption of zero error
% cross-covariance.
E1 = Q11 - (Q12 * Q13) / Q23;
E2 = Q22 - (Q12 * Q23) / Q13;
E3 = Q33 - (Q13 * Q23) / Q12;

Evec = [E1; E2; E3];

%% ---------- Strict validity condition ----------
% If any estimated error variance is negative (or invalid), the entire
% TC solution is considered physically inconsistent.
if any(Evec < 0) || any(isnan(Evec)) || any(isinf(Evec))
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

% Error covariance matrix is diagonal by construction.
EeeT_est = diag(Evec);

%% ---------- Signal variance estimation ----------
% The latent signal variance can be estimated from any cyclic relation.
Ey2_est = (Q12 * Q13) / Q23;

% If Ey2 becomes invalid, the full solution is invalid.
if isnan(Ey2_est) || isinf(Ey2_est) || Ey2_est <= 0
    EeeT_est  = NaN(p);
    theta_est = NaN(p,1);
    rho2_est  = NaN(p,1);
    return
end

%% ---------- Scaling factor estimation ----------
% theta is only identifiable up to a multiplicative constant unless a
% reference normalization is adopted. Here the first dataset is used as
% the reference scale, giving theta_1 = 1.
theta_est = [1;
             Q23 / Q13;
             Q23 / Q12];

%% ---------- Squared correlation estimation ----------
rho2_est = 1 - Evec / Ey2_est;

% Numerical cleanup:
% tiny roundoff violations outside [0,1] are clipped back.
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;

end