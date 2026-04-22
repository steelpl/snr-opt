function [N_est,a_est] = SNRest(ExxT, Ey2)
% ============================================================
% SNR-based decomposition (SNRest)
% Rank-1 signal extraction from normalized covariance
%
% INPUT
%   ExxT : covariance matrix of observations (p x p)
%   Ey2  : signal power (scalar)
%
% OUTPUT
%   N_est : estimated noise-to-signal ratio matrix (p x p)
%   a_est : estimated scaling factor vector (p x 1)
%
% METHOD
%   The method operates on the normalized covariance:
%
%       C = ExxT / Ey2
%
%   and estimates a rank-1 signal structure:
%
%       C ≈ a a^T + N
%
%   by minimizing the off-diagonal mismatch of (C - a a^T)
%   using projected subgradient descent.
%
%   The constraint
%
%       a_i^2 ≤ C_ii
%
%   ensures non-negative noise variances:
%
%       diag(N_est) = C_ii - a_i^2 ≥ 0
%
% RELATION TO NC FRAMEWORK
%   SNRest is the normalized counterpart of the NC formulation:
%
%       Q = θ θ^T + Σ_e
%
%   with the relation:
%
%       θ = a * sqrt(Ey2)
%
%   Thus, SNRest can be interpreted as a special case of NC
%   applied to normalized covariance space.
%
% NOTES
%   - The solution is not unique in scale; it is tied to Ey2.
%   - The algorithm is robust to moderate error cross-correlation.
%   - Projection guarantees physically feasible diagonal estimates.
%
% REFERENCE
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Trans Geosci Remote Sens.
%
%   If you use this method, please cite the above reference.
% ============================================================

%% ---------- Input check ----------
if nargin < 2
    error('SNRest requires two inputs: ExxT and Ey2.');
end

[p1,p2] = size(ExxT);
if p1 ~= p2
    error('ExxT must be a square matrix.');
end

if ~isscalar(Ey2) || Ey2 <= 0
    error('Ey2 must be a positive scalar.');
end

%% ---------- Setup ----------
p = p1;

% Symmetrize for numerical stability
Q = (ExxT + ExxT') / 2;

% Normalize covariance
C = Q / Ey2;

diagC = diag(C);
if any(diagC <= 0)
    error('Normalized covariance has non-positive diagonal entries.');
end

sqrtDiagC = sqrt(diagC);

%% ---------- Parameters ----------
beta   = 0.5 * min(diagC);   % isotropic noise approximation
lambda = 0.01;               % step size
iters  = 2000;
tol    = 1e-12;

%% ---------- Initialization ----------
% Leading eigenpair gives best rank-1 approximation under L2 sense
opts.tol = 1e-6;
opts.maxit = 500;

[V,D] = eigs(C - beta*eye(p), 1, 'largestreal', opts);

a_est = V * sqrt(max(D,0));

% Enforce consistent sign (stabilizes iteration)
a_est = sign(a_est) .* abs(a_est);

% Scale step size relative to magnitude
lambda = lambda / max(norm(a_est),1e-8);

%% ---------- Iterative optimization ----------
for i = 1:iters
    
    % Rank-1 reconstruction
    A = a_est * a_est';
    
    % Subgradient of L1 off-diagonal mismatch
    S = sign(A - C);
    S(1:p+1:end) = 0;   % remove diagonal contribution
    
    grad = S * a_est;
    
    % Gradient descent
    a_est = a_est - lambda * grad;
    
    % Projection: enforce a_i^2 <= C_ii
    a_est = sign(a_est) .* min(abs(a_est), sqrtDiagC);
    
end

%% ---------- Final outputs ----------
A = a_est * a_est';

% Noise-to-signal ratio matrix
N_est = C - A;

% Symmetry cleanup
N_est = (N_est + N_est') / 2;

% Numerical cleanup for diagonal
d = diag(N_est);
d(d < 0 & d > -tol) = 0;
N_est(1:p+1:end) = d;

end