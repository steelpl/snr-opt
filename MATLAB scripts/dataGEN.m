function [y,e] = dataGEN(n,p,ecc,SNRdB)
% ============================================================
% Synthetic data generator for signal and error components
%
% INPUT
%   n     : data length (scalar)
%   p     : number of datasets (scalar)
%   ecc   : error cross-correlation (scalar, [0,1])
%   SNRdB : signal-to-noise ratio in dB (scalar)
%
% OUTPUT
%   y : signal (n x 1), zero mean
%   e : error  (n x p), zero mean
%
% METHOD
%   The function generates synthetic signal and error samples such that
%   the joint covariance of [y, e] matches a prescribed target structure:
%
%       cov([y, e]) = [ Ey2     0
%                        0    EeeT ]
%
%   where EeeT is a user-controlled error covariance matrix with
%   diagonal variances and off-diagonal correlations governed by ecc.
%
% NOTES
%   - The generated y and e are uncorrelated by construction.
%   - Positive semi-definiteness is enforced numerically for both EeeT
%     and the joint covariance matrix.
%   - The final samples are obtained by whitening Gaussian noise and
%     recoloring it to match the target covariance.
%
% REFERENCE
% For more details, see:
%
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Trans Geosci Remote Sens., 60, 1–15.
%
% If you use the methods presented in the paper and/or this example, 
% please cite this paper where appropriate.
% ============================================================

%% ---------- Input check ----------
if nargin < 4
    error('dataGEN requires four inputs: n, p, ecc, and SNRdB.');
end

if ~isscalar(n) || n <= 1 || floor(n) ~= n
    error('n must be an integer greater than 1.');
end

if ~isscalar(p) || p < 1 || floor(p) ~= p
    error('p must be a positive integer.');
end

if ~isscalar(SNRdB)
    error('SNRdB must be a scalar.');
end

%% ---------- Setup ----------
SNR = 10^(SNRdB/10);
ecc = max(0, min(1, ecc));   % keep ecc within [0,1]
tol = 1e-10;

b = nchoosek(1:p,2); % combination pairs e.g. for p=3, (1,2), (1,3), (2,3)

%% ---------- Error covariance construction ----------
% Random error variances for each dataset.
var_e = rand(p,1);

% Initialize EeeT with diagonal variances.
EeeT = diag(var_e);

% Set pairwise error covariance using a common correlation coefficient ecc.
for i = 1:size(b,1)
    ii = b(i,1);
    jj = b(i,2);
    
    val = sqrt(var_e(ii)*var_e(jj)) * ecc;
    EeeT(ii,jj) = val;
    EeeT(jj,ii) = val;
end

% Symmetry cleanup and numerical PSD enforcement.
EeeT = (EeeT + EeeT') / 2;
[V,D] = eig(EeeT);
D = max(D, tol*eye(p));
EeeT = V * D * V';

%% ---------- Target joint covariance ----------
% Signal power is chosen to match the prescribed SNR relative to the
% average error variance.
Ey2 = mean(diag(EeeT) * SNR);

% Joint covariance matrix of [y, e].
m = zeros(p+1);
m(1,1) = Ey2;
m(2:end,2:end) = EeeT;

% Symmetry cleanup and numerical PSD enforcement for the joint matrix.
m = (m + m') / 2;
[V,D] = eig(m);
D = max(D, tol*eye(p+1));
m = V * D * V';

%% ---------- Generate Gaussian samples ----------
% Start from standard Gaussian samples.
ye = randn(n, p+1); % [y, e]

% Remove sample means column-wise so that the generated samples are
% centered before covariance shaping.
ye = ye - mean(ye,1);

%% ---------- Whitening ----------
% Compute sample covariance of the raw Gaussian samples.
C = cov(ye);
C = (C + C') / 2;

% Enforce numerical PSD for stable Cholesky factorization.
[Vc,Dc] = eig(C);
Dc = max(Dc, tol*eye(p+1));
C = Vc * Dc * Vc';

% Whiten the samples so that the covariance becomes approximately identity.
ye = ye / chol(C);

%% ---------- Recoloring ----------
% Impose the target covariance structure m.
ye = ye * chol(m);

%% ---------- Output ----------
y = ye(:,1);
e = ye(:,2:end);

end
