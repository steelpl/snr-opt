function u = maxR(theta, Q)
% ============================================================
% maxR: Correlation-maximizing data fusion weights
%
% INPUT
%   theta : signal scaling vector (p x 1)
%           (typically obtained from NC or SNR-based estimation)
%   Q     : covariance matrix of observations (p x p)
%
% OUTPUT
%   u     : optimal merging weights (p x 1)
%
% METHOD
%   The weights are obtained by maximizing the squared Pearson
%   correlation between the merged estimate (u'x) and the latent signal.
%
%   This leads to the generalized Rayleigh quotient:
%
%       max_u   (u' θθ' u) / (u' Q u)
%
%   whose solution satisfies:
%
%       Q u = θ
%
%   Thus, the optimal weights are given by:
%
%       u ∝ Q^{-1} θ
%
% NOTES
%   - The solution is defined up to a multiplicative constant.
%   - Normalization (sum(u)=1) is applied for interpretability,
%     but is not required for optimality.
%   - The pseudoinverse handles singular observation covariances.
%
% REFERENCE
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Transactions on Geoscience and Remote Sensing, 60, 1–15.
%   https://doi.org/10.1109/TGRS.2021.3107028
%
% Also related to classical Rayleigh quotient optimization:
% Horn, R. A., & Johnson, C. R. (2013).
% Matrix Analysis (2nd ed.). Cambridge University Press.
% ============================================================

%% ---------- Input check ----------
if nargin < 2
    error('maxR requires two inputs: theta and Q.');
end

if ~isvector(theta)
    error('theta must be a vector.');
end

theta = theta(:); % enforce column vector

[p1,p2] = size(Q);
if p1 ~= p2
    error('Q must be a square matrix.');
end

if length(theta) ~= p1
    error('Dimensions of theta and Q do not match.');
end

%% ---------- Setup ----------
% Symmetrize Q to remove small numerical asymmetry.
Q = 0.5 * double(Q) + 0.5 * double(Q');
theta = double(theta);

%% ---------- Solve optimal weights ----------
% Solve Q u = theta (Rayleigh quotient solution)
% Using backslash ensures numerical stability over explicit inverse.
u = pinv(Q) * theta;

%% ---------- Normalization ----------
% Normalize weights for interpretability (sum to 1).
% This does not affect optimality.
s = sum(u);

if abs(s) > 1e-12 * norm(u,1)
    u = u / s;
end

end
