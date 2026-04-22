function [EeeT_est, theta_est, rho2_est] = NC(covx)
% ============================================================
% N-Tuple Collocation (NC)
% Constrained rank-1 approximation of covariance matrix
%
% INPUT
%   covx : covariance matrix (p x p)
%
% OUTPUT
%   EeeT_est  : estimated error covariance matrix (p x p)
%   theta_est : estimated scaling vector (p x 1)
%   rho2_est  : squared correlation with latent signal (p x 1)
%
% METHOD
%   The method approximates
%       Q = theta*theta' + EeeT
%   by minimizing the off-diagonal residuals of (Q - theta*theta')
%   under the constraint theta_i^2 <= Q_ii, which guarantees
%   non-negative diagonal error variances.
%
% NOTES
%   - The objective is based on an L1-type off-diagonal mismatch.
%   - The projection step enforces physical feasibility.
%   - For p = 3, the solution is expected to be close to classical TC.
% ============================================================

%% ---------- Input check ----------
if nargin < 1
    error('NC requires one input argument: covx');
end

if ~ismatrix(covx) || size(covx,1) ~= size(covx,2)
    error('Input covx must be a square matrix.');
end

%% ---------- Setup ----------
p = size(covx,1);

% Symmetrize Q to remove tiny numerical asymmetry in the input.
Q = (covx + covx') / 2;

diagQ = diag(Q);

% If any diagonal entry is negative, the covariance matrix is invalid.
if any(diagQ < 0)
    error('Input covariance matrix has negative diagonal entries.');
end

% Precompute square roots of diagonal entries for projection.
sqrtDiagQ = sqrt(max(diagQ,0));

%% ---------- Parameters ----------
beta = 0.5 * min(diagQ);   % isotropic error approximation for initialization
eta   = 0.01;              % nominal step size
iters = 2000;              % number of iterations
tol   = 1e-12;             % numerical tolerance for final cleanup

% Mask that keeps only off-diagonal elements in the objective.
mask = ones(p) - eye(p);

%% ---------- Spectral Initialization ----------
% Initial theta is obtained from the leading eigenpair of (Q - beta I),
% which approximates a rank-1 signal component under isotropic error.
[vecs, vals] = eig(Q - beta*eye(p));
[lambda1, idx] = max(diag(vals));
v1 = vecs(:,idx);

theta_est = sqrt(max(lambda1,0)) * v1;

% Sign convention: make theta consistently nonnegative where possible.
% This does not change theta*theta', but stabilizes the iterative updates.
theta_est = sign(theta_est) .* abs(theta_est);

% Scale the step size by the norm of the initial solution so that
% the effective update magnitude is less sensitive to the data scale.
eta = eta / max(norm(theta_est),1e-8);

%% ---------- Iterative Optimization ----------
for k = 1:iters

    % Current rank-1 signal covariance estimate.
    theta_outer = theta_est * theta_est';

    % Subgradient of the L1 off-diagonal mismatch:
    %   sum_{i ~= j} |Q_ij - theta_i theta_j|
    grad = (mask .* sign(theta_outer - Q)) * theta_est;

    % Gradient descent update.
    theta_est = theta_est - eta * grad;

    % Projection step:
    % enforce theta_i^2 <= Q_ii so that diag(EeeT_est) = Q_ii - theta_i^2 >= 0.
    theta_est = sign(theta_est) .* min(abs(theta_est), sqrtDiagQ);

end

%% ---------- Final Outputs ----------
theta_outer = theta_est * theta_est';

% Estimated error covariance matrix.
EeeT_est = Q - theta_outer;

% Symmetry cleanup to suppress tiny asymmetry from floating-point arithmetic.
EeeT_est = (EeeT_est + EeeT_est') / 2;

% Numerical cleanup for diagonal entries:
% tiny negative values caused only by roundoff are set to zero.
d = diag(EeeT_est);
d(d < 0 & d > -tol) = 0;
EeeT_est(1:p+1:end) = d;

% Squared correlation with the latent signal.
rho2_est = (theta_est.^2) ./ max(diagQ, tol);

% Optional numerical cleanup:
% rho2 should theoretically lie in [0,1] when the model is valid.
rho2_est(rho2_est < 0 & rho2_est > -tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + tol) = 1;

end