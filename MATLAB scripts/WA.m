function u = WA(EeeT)
% ============================================================
% WA: Error-covariance-weighted average
%
% INPUT
%   EeeT : error covariance matrix (p x p)
%
% OUTPUT
%   u    : merging weight (p x 1)
%
% METHOD
%   The weighted-average solution is obtained by minimizing the variance
%   of the merged estimate under the unbiasedness constraint:
%
%       min_u   u' EeeT u
%       s.t.    eta' u = 1
%
%   where eta is a p x 1 vector of ones.
%
%   The closed-form solution is:
%
%       u = EeeT^(-1) eta / (eta' EeeT^(-1) eta)
%
%   In implementation, MATLAB backslash is used instead of an explicit
%   matrix inverse for numerical stability.
%
% NOTES
%   - This is the classical variance-minimizing weighted average.
%   - The method uses only the error covariance matrix and does not
%     explicitly account for scaling factors.
%   - For this reason, WA is most directly comparable to other methods
%     when the input datasets are on a common scale.
%
% REFERENCE
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Trans Geosci Remote Sens.
%
%   If you use this method, please cite the above reference.
% ============================================================

%% ---------- Input check ----------
if nargin < 1
    error('WA requires one input: EeeT.');
end

[p1,p2] = size(EeeT);
if p1 ~= p2
    error('EeeT must be a square matrix.');
end

%% ---------- Setup ----------
% Symmetrize EeeT to suppress tiny numerical asymmetry.
EeeT = (EeeT + EeeT') / 2;

eta = ones(p1,1);

%% ---------- Solve optimal weights ----------
% Compute EeeT^{-1} eta using backslash.
tmp = EeeT \ eta;

% Normalize so that sum(u) = 1.
u = tmp / (eta' * tmp);

end