function u = WA(EeeT, covarianceScale)
%WA Minimum-error-variance weights subject to sum(u) = 1.
%   EeeT must be a positive semidefinite error covariance matrix. Products
%   should already be on a common scale. For singular EeeT, use a
%   zero-error combination if the constraint can be satisfied in its
%   nullspace. Otherwise use the minimum-norm constrained solution.
%   WA(E, covarianceScale) uses the observation covariance scale in the
%   roundoff threshold for E formed by subtraction in NC. For an NC fit,
%   pass fit.info.covarianceScale.
%
%   Kim et al. (2022). IEEE Transactions on Geoscience and Remote Sensing, 60.
%   https://doi.org/10.1109/TGRS.2021.3107028

validateattributes(EeeT, {'numeric'}, {'real','finite','square','nonempty'}, ...
    'WA', 'EeeT', 1);
EeeT = double(EeeT);
[V, d] = eig(0.5 * EeeT + 0.5 * EeeT', 'vector');
% Treat roundoff near the PSD boundary, including NC's interior margin, as zero.
eigenScale = max(abs(d));
if nargin >= 2
    validateattributes(covarianceScale, {'numeric'}, {'real','finite','scalar','nonnegative'}, ...
        'WA', 'covarianceScale', 2);
    eigenScale = max(eigenScale,double(covarianceScale));
end
tol = 64 * numel(d) * eps * eigenScale;
if any(d < -tol)
    error('WA:IndefiniteCovariance', 'Error covariance must be positive semidefinite.');
end

onesVector = ones(numel(d),1);
positive = d > tol;
nullComponent = V(:,~positive)' * onesVector;
if norm(nullComponent) > 64 * numel(d) * eps * norm(onesVector)
    direction = V(:,~positive) * nullComponent;
else
    % Rescale eigenvalues before taking reciprocals to avoid an overflow.
    direction = V(:,positive) * ((V(:,positive)' * onesVector) .* ...
        (min(d(positive)) ./ d(positive)));
end
u = direction / sum(direction);
end
