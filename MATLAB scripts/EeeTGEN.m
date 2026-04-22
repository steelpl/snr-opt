function EeeT = EeeTGEN(p,ecc,var_e)
% ============================================================
% Synthetic error covariance matrix generator
%
% INPUT
%   p     : number of datasets (scalar)
%   ecc   : error cross-correlation (scalar, [0,1])
%   var_e : error variances (p x 1, optional)
%
% OUTPUT
%   EeeT : error covariance matrix (p x p)
%
% METHOD
%   The function constructs a covariance matrix whose diagonal entries are
%   given by var_e and whose off-diagonal entries are controlled by a
%   common correlation coefficient ecc:
%
%       EeeT(i,j) = sqrt(var_e(i)*var_e(j)) * ecc,   i ~= j
%
%   If var_e is not provided, random positive variances are generated.
%
% NOTES
%   - The resulting matrix is symmetrized explicitly.
%   - Numerical positive semi-definiteness is enforced at the end to avoid
%     issues caused by floating-point roundoff.
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
if nargin < 2
    error('EeeTGEN requires at least two inputs: p and ecc.');
end

if ~isscalar(p) || p < 1 || floor(p) ~= p
    error('p must be a positive integer.');
end

% Keep ecc within the intended range.
ecc = max(0, min(1, ecc));

if nargin < 3 || isempty(var_e)
    var_e = rand(p,1);
else
    if ~isvector(var_e) || length(var_e) ~= p
        error('var_e must be a vector of length p.');
    end
    var_e = var_e(:);
    
    if any(var_e < 0)
        error('var_e must contain non-negative entries.');
    end
end

tol = 1e-10;

%% ---------- Covariance construction ----------
% Combination pairs, e.g. for p = 3:
% (1,2), (1,3), (2,3)
b = nchoosek(1:p,2);

% Initialize with diagonal error variances.
EeeT = diag(var_e);

% Fill off-diagonal entries using the prescribed common
% error cross-correlation coefficient.
for i = 1:size(b,1)
    ii = b(i,1);
    jj = b(i,2);

    val = sqrt(var_e(ii) * var_e(jj)) * ecc;

    EeeT(ii,jj) = val;
    EeeT(jj,ii) = val;
end

%% ---------- Numerical cleanup ----------
% Enforce symmetry explicitly to suppress tiny asymmetry caused by
% floating-point arithmetic.
EeeT = (EeeT + EeeT') / 2;

% Enforce numerical positive semi-definiteness.
[V,D] = eig(EeeT);
D = max(D, tol*eye(p));
EeeT = V * D * V';

end
