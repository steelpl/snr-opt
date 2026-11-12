function [N_est, a_est, snr_est, fit] = SNRest(ExxT, Ey2, opts)
%SNREST Estimate SNR and merging parameters from one NC covariance fit.
%   [N,A] = SNRest(Q,Ey2) preserves the normalized decomposition interface.
%   Q is the covariance of collocated, centered products. Ey2 is the
%   positive target signal variance. Omit it or pass [] for a unit-variance
%   latent target. OPTS is forwarded to NC.
%
%   [N,A,SNR,FIT] also returns the linear SNR per product and NC outputs.
%   FIT contains errorCov, theta, rho2, covariance (used by NC), and info.
%   N = errorCov/Ey2, A = theta/sqrt(Ey2), SNR = theta.^2./diag(errorCov).
%   The loading with largest magnitude is oriented positive by convention.
%   To align to a known opposite direction, flip A and FIT.theta together.
%   Zero error with positive signal gives Inf. Zero total variance gives NaN.
%
%   Pham, H. T., & Kim, S. (2027). Journal of Hydrology, 682, 136749.
%   https://doi.org/10.1016/j.jhydrol.2026.136749
%   Kim et al. (2022). IEEE Transactions on Geoscience and Remote Sensing, 60.
%   https://doi.org/10.1109/TGRS.2021.3107028

if nargin < 2 || isempty(Ey2)
    Ey2 = 1;
end
validateattributes(Ey2, {'numeric'}, {'real','scalar','finite','positive'}, ...
    'SNRest', 'Ey2', 2);
Ey2 = double(Ey2);
if nargin < 3
    opts = [];
end

[errorCov, theta, rho2, info] = NC(ExxT, opts);
[~, anchor] = max(abs(theta));
if theta(anchor) < 0
    theta = -theta;
end

N_est = errorCov / Ey2;
a_est = theta / sqrt(Ey2);
signalVariance = theta.^2;
noiseVariance = max(diag(errorCov), 0);
snr_est = signalVariance ./ noiseVariance;
fit = struct('errorCov', errorCov, 'theta', theta, 'rho2', rho2, ...
    'covariance', errorCov + theta*theta', 'info', info);
end
