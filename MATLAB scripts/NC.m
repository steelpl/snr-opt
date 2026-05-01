function [EeeT_est, theta_est, rho2_est] = NC(covx)
% N-Tuple Collocation (NC)
% Constrained rank-1 decomposition of a covariance matrix
%
% INPUT
%   covx : p-by-p covariance matrix
%
% OUTPUT
%   EeeT_est  : estimated error covariance (p-by-p), symmetric
%   theta_est : estimated latent-scaling vector (p-by1)
%   rho2_est  : squared correlation with latent signal (p-by1)
%
% PURPOSE
%   Find theta and EeeT satisfying
%       Q ≈ theta*theta' + EeeT
%   with the constraints theta_i^2 <= Q_ii so EeeT has nonnegative diagonals.
%   The fit minimizes off-diagonal mismatch (L1-like) between Q and theta*theta'.
%
% METHOD SUMMARY
%   - Symmetrize input Q and initialize theta from largest eigenpair of Q - beta*I.
%   - Use iterative subgradient descent on off-diagonal absolute deviations:
%         grad ≈ (sign(theta*theta' - Q) with zeroed diagonal) * theta
%   - After each step project each theta_i to satisfy |theta_i| <= sqrt(Q_ii).
%   - Stop on small relative change or max iterations.
%   - Return symmetric residual EeeT_est = Q - theta*theta' and elementwise
%     bounded rho2_est = theta_i^2 ./ Q_ii.
%
% IMPLEMENTATION NOTES
%   - Small tolerances guard against tiny negative numerics on diagonals.
%   - Masking keeps diagonal entries from influencing the off-diagonal objective.
% ============================================================

% Validate input presence and shape
if nargin < 1
    error('NC requires one input argument: covx');
end
if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end
p = size(covx,1);

% Symmetrize input and check diagonal feasibility
Q = 0.5 * (covx + covx');
diagQ = diag(Q);
if any(diagQ < 0)
    error('Input covariance matrix has negative diagonal entries.');
end
sqrtDiagQ = sqrt(max(diagQ, 0));

% Parameters (adjusted defaults)
beta     = 0.5 * min(diagQ);
eta      = 0.01;
maxIters = 2000;
tol      = 1e-12;
relTol   = 1e-8;    % relative change tolerance for early stopping

% Mask for off-diagonals (logical and sparse to speed heavy ops when p large)
mask = sparse(ones(p) - eye(p));

% Spectral initialization
[vecs, vals] = eig(Q - beta * eye(p));
[~, idx]     = max(diag(vals));
v1           = vecs(:, idx);
theta_est    = sqrt(max(vals(idx, idx), 0)) * v1;
theta_est    = sign(theta_est) .* abs(theta_est);    % keep sign convention

% Scale step size relative to initialization magnitude
eta = eta / max(norm(theta_est), 1e-8);

% Iterative optimization (subgradient descent with projection)
theta_old = theta_est;
for k = 1:maxIters
    % Outer product (full) needed for sign comparisons
    theta_outer = theta_est * theta_est';
    
    % Off-diagonal sign subgradient (choose 0 for exact zeros)
    S = sign(theta_outer - Q);  % p-by-p
    S = S .* mask;              % zero diagonal
    
    % Compute gradient and take a step
    grad = S * theta_est;       % p-by-1
    theta_est = theta_est - eta * grad;
    
    % Projection: enforce theta_i^2 <= Q_ii (non-negative error variances)
    abs_theta = abs(theta_est);
    abs_theta = min(abs_theta, sqrtDiagQ);
    theta_est = sign(theta_est) .* abs_theta;
    
    % Early stopping based on relative change
    if norm(theta_est - theta_old) <= relTol * max(1, norm(theta_old))
        break;
    end
    theta_old = theta_est;
end

% Form residual estimate and sanitize numerical noise
theta_outer = theta_est * theta_est';
EeeT_est    = Q - theta_outer;
EeeT_est    = 0.5 * (EeeT_est + EeeT_est');  % ensure symmetry

d = diag(EeeT_est);
d(d < 0 & d > -tol) = 0;
EeeT_est(1:p+1:end) = d;

rho2_est = (theta_est .^ 2) ./ max(diagQ, tol);
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;
end