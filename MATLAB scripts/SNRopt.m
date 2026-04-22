function u = SNRopt(N,a)
% ============================================================
% SNRopt: Signal-to-noise-ratio-optimal merging weights
%
% INPUT
%   N : noise-to-signal ratio matrix (p x p)
%   a : scaling factor vector (p x 1)
%
% OUTPUT
%   u : optimal merging weights (p x 1)
%
% METHOD
%   The SNR-optimal weights are obtained from the linear system:
%
%       (N + a a^T) u = a
%
%   where
%       N      : noise-to-signal ratio matrix
%       a a^T  : rank-1 signal covariance in normalized space
%
%   Thus, the solution is:
%
%       u = (N + a a^T)^(-1) a
%
%   In implementation, the system is solved using MATLAB backslash
%   rather than an explicit matrix inverse for numerical stability.
%
% NOTES
%   - The formulation assumes that the signal and noise decomposition
%     is expressed in normalized covariance space:
%
%         C = a a^T + N
%
%   - This method is directly linked to the SNRest framework and can
%     be interpreted as the normalized-space counterpart of maxR.
%
% REFERENCE
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%   Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%   IEEE Trans Geosci Remote Sens., 60, 1–15.
%
%   If you use this method, please cite the above reference.
% ============================================================

%% ---------- Input check ----------
if nargin < 2
    error('SNRopt requires two inputs: N and a.');
end

[p1,p2] = size(N);
if p1 ~= p2
    error('N must be a square matrix.');
end

if ~isvector(a)
    error('a must be a vector.');
end

a = a(:); % enforce column vector

if length(a) ~= p1
    error('Dimensions of N and a do not match.');
end

%% ---------- Setup ----------
% Symmetrize N to suppress tiny numerical asymmetry.
N = (N + N') / 2;

%% ---------- Solve optimal weights ----------
% Solve (N + a a^T) u = a
% Using backslash is more stable than explicit inversion.
u = (N + a*a') \ a;

end
