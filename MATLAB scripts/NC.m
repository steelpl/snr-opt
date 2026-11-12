function [EeeT_est, theta_est, rho2_est, info] = NC(covx, opts)
% ============================================================
% N-Tuple Collocation (NC)
% Algorithm from https://github.com/steelpl/NC-public (be461e3, MIT).
% Constrained rank-1 decomposition of a covariance matrix
%
% INPUT
%   covx : p-by-p covariance matrix
%   opts : optional structure controlling numerical optimization
%          opts.nStarts       : number of deterministic spectral starts
%          opts.nRandomStarts : number of additional random starts
%          opts.randomSeed    : random seed for random starts
%          opts.maxIters      : maximum iterations per start
%          opts.patience      : early-stopping patience
%          opts.enforcePSD    : enforce PSD error covariance, default true
%                              false is a raw diagnostic comparison only
%          opts.repairInputPSD: repair indefinite input covariance, default true
%          opts.verbose       : print optimization summary
%
% OUTPUT
%   EeeT_est  : estimated error covariance (p-by-p), symmetric
%   theta_est : estimated latent-scaling vector (p-by-1)
%   rho2_est  : squared correlation with latent signal (p-by-1)
%   info      : diagnostic structure with objective and convergence details
%
% PURPOSE
%   Find theta and EeeT satisfying
%       Q ≈ theta*theta' + EeeT
%   First fit with theta_i^2 <= Q_ii, then (by default) project theta into
%   range(Q) and shrink it radially until theta' * pinv(Q) * theta <= 1.
%   This makes Q - theta*theta' PSD while retaining the decomposition.
%   This two-stage estimator is NOT an exact solution of a joint PSD-
%   constrained L1 optimization. The loss change is reported explicitly.
%
%   The fit minimizes the off-diagonal L1 mismatch between Q and
%   theta*theta':
%
%       sum_{i ~= j} |Q_ij - theta_i theta_j|
%
% METHOD SUMMARY
%   - Symmetrize input covariance matrix Q.
%   - Work in globally normalized covariance units for scale invariance.
%   - Initialize from Q - beta*I, with beta <= 0.95*lambda_max(Q).
%   - Apply projected subgradient descent with monotone backtracking to
%     reduce the off-diagonal L1 mismatch.
%   - After each step, project theta_i to satisfy |theta_i| <= sqrt(Q_ii).
%   - Stop adaptively when the off-diagonal L1 objective no longer improves.
%   - Return the best theta found across all starts.
%   - By default, apply the PSD feasibility map to the final theta.
%
% NOTES
%   - maxIters is only a safety cap.
%   - Early stopping is controlled by objective improvement, not iteration count.
%   - Small tolerances guard against tiny negative numerical values.
%   - A small objective does not identify the unique physical error matrix.
%   - p = 2 is computationally supported but structurally underdetermined;
%     p >= 3 is recommended and still requires assumptions/regularization.
%   - Cov(true signal,error_i) = 0 is needed to interpret the decomposition.
%
% REFERENCES
%   Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).
%     Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
%     IEEE Transactions on Geoscience and Remote Sensing, 60, 1–15.
%     https://doi.org/10.1109/TGRS.2021.3107028
%
%   Pham, H. T., & Kim, S. (2027).
%     N-Tuple Collocation for Reference-Free Uncertainty Estimation from
%     Multi-Source Hydrological Data: Beyond Triple Collocation and
%     Three-Cornered Hat.
%     Journal of Hydrology, 682, 136749.
%     https://doi.org/10.1016/j.jhydrol.2026.136749
%
% If you use the methods presented in the paper and/or this example,
% please cite this paper where appropriate.
% ============================================================

%% ---------- Input validation ----------
if nargin < 1
    error('NC requires one input argument: covx');
end

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ~isnumeric(covx) || ~isreal(covx) || ~ismatrix(covx) || ...
        size(covx,1) ~= size(covx,2) || size(covx,1) < 2
    error('NC:InvalidCovariance', 'Input covx must be a real square matrix with p >= 2.');
end

if any(~isfinite(covx(:)))
    error('Input covx must contain only finite values.');
end

p = size(covx,1);
if ~isstruct(opts) || ~isscalar(opts)
    error('NC:InvalidOptions', 'opts must be a scalar structure.');
end

%% ---------- Parameters ----------
params.betaFactor    = getOpt(opts, 'betaFactor', 0.5);
% Half the former default step because the exact symmetric gradient is 2*S*theta.
params.eta0          = getOpt(opts, 'eta0', 0.005);
params.maxIters      = getOpt(opts, 'maxIters', 2000);
params.minIters      = getOpt(opts, 'minIters', 20);
params.patience      = getOpt(opts, 'patience', 30);
params.tol           = getOpt(opts, 'tol', 1e-12);
params.objTol        = getOpt(opts, 'objTol', 1e-10);
params.etaDecay      = getOpt(opts, 'etaDecay', 0.7);
params.etaMin        = getOpt(opts, 'etaMin', 1e-10);
params.lineSearchMax = getOpt(opts, 'lineSearchMax', 12);
params.nStarts       = min(p, max(1, getOpt(opts, 'nStarts', min(4, p))));
params.nRandomStarts = max(0, getOpt(opts, 'nRandomStarts', 0));
params.randomSeed    = getOpt(opts, 'randomSeed', 1);
params.enforcePSD    = getOpt(opts, 'enforcePSD', true);
params.repairInputPSD = getOpt(opts, 'repairInputPSD', true);
params.psdTol        = getOpt(opts, 'psdTol', 1e-10);
params.verbose       = getOpt(opts, 'verbose', false);
validate_parameters(params);

%% ---------- Symmetrize and check covariance matrix ----------
Qoriginal = 0.5 * double(covx) + 0.5 * double(covx');
covarianceScale = max(abs(Qoriginal(:)));
if covarianceScale == 0
    covarianceScale = 1; % exact zero input; no unit-dependent positive floor
end
Q = Qoriginal / covarianceScale;
[Q, inputPSDInfo] = repair_psd_matrix(Q, params.psdTol, params.repairInputPSD, ...
    'Input covariance matrix');
diagQ = diag(Q);

if any(diagQ < -params.tol)
    error('Input covariance matrix has negative diagonal entries.');
end

diagQ(diagQ < 0 & diagQ > -params.tol) = 0;
Q(1:p+1:end) = diagQ;
sqrtDiagQ = sqrt(max(diagQ,0)); % projection bound for theta

requestedBeta = params.betaFactor * min(diagQ);
maxInputEigenvalue = max(eig(Q));
beta = min(max(requestedBeta, 0), 0.95 * max(maxInputEigenvalue, 0));

%% ---------- Off-diagonal mask ----------
offMask = true(p);
offMask(1:p+1:end) = false;     % exclude diagonal from L1 objective

%% ---------- Multi-start projected subgradient descent ----------
thetaStarts = build_theta_starts(Q,beta,sqrtDiagQ,params);

bestGlobalObj = Inf;
bestGlobalTheta = NaN(p,1);
startInfo = repmat(struct( ...
    'initialObjective', NaN, ...
    'finalObjective', NaN, ...
    'iterations', 0, ...
    'converged', false), size(thetaStarts,2), 1);

for iStart = 1:size(thetaStarts,2)

    [theta_i, obj_i, runInfo] = run_projected_l1( ...
        Q, thetaStarts(:,iStart), sqrtDiagQ, offMask, params);

    startInfo(iStart) = runInfo;

    if obj_i < bestGlobalObj
        bestGlobalObj = obj_i;
        bestGlobalTheta = theta_i;
        bestStart = iStart;
    end
end

%% ---------- Use best solution ----------
theta_est = bestGlobalTheta;

psdInfo = struct('applied', false, 'thetaScale', 1, ...
    'rangeProjectionRelativeChange', 0, 'thetaRelativeChange', 0, ...
    'mahalanobisSquaredBeforeShrink', NaN, ...
    'minEigenvalueBefore', NaN, 'minEigenvalueAfter', NaN);
if params.enforcePSD
    [theta_est, psdInfo] = enforce_residual_psd(Q, theta_est, params.psdTol);
end

theta_outer = theta_est * theta_est';
EeeT_est = Q - theta_outer;
EeeT_est = 0.5 * (EeeT_est + EeeT_est');   % enforce symmetry
% Do not independently clip E: doing so would break Q = theta*theta' + E.
finalPSDInfo = psd_matrix_info(EeeT_est);

%% ---------- Squared correlation ----------
rho2_est = NaN(p,1); % correlation is undefined for a zero-variance dataset
positiveVariance = diagQ > 0;
rho2_est(positiveVariance) = theta_est(positiveVariance).^2 ./ diagQ(positiveVariance);

rho2_est(rho2_est < 0 & rho2_est > -params.tol) = 0;
rho2_est(rho2_est > 1 & rho2_est < 1 + params.tol) = 1;

finalObjective = offdiag_L1_objective(Q, theta_est, offMask);

%% ---------- Diagnostics ----------
info = struct();
info.objective = finalObjective * covarianceScale;
info.objectiveBeforePSD = bestGlobalObj * covarianceScale;
info.objectiveChangePSD = (finalObjective - bestGlobalObj) * covarianceScale;
info.objectiveChangePSDRelativeToQ = (finalObjective - bestGlobalObj) / ...
    max(norm(Q,'fro'), realmin);
info.bestStart = bestStart;
info.nStarts = size(thetaStarts,2);
info.startInfo = startInfo;
for iStart = 1:numel(info.startInfo)
    info.startInfo(iStart).initialObjective = startInfo(iStart).initialObjective * covarianceScale;
    info.startInfo(iStart).finalObjective = startInfo(iStart).finalObjective * covarianceScale;
end
info.minErrorVariance = min(diag(EeeT_est)) * covarianceScale;
info.minErrorCovEigenvalue = finalPSDInfo.minEigenvalueAfter * covarianceScale;
info.maxRho2 = max(rho2_est);
info.feasibleDiagonal = all(diag(EeeT_est) >= -params.tol);
info.feasiblePSD = finalPSDInfo.minEigenvalueAfter >= -params.psdTol * max(norm(Q,'fro'),realmin);
info.inputMinEigenvalue = inputPSDInfo.minEigenvalueBefore * covarianceScale;
info.inputPSDRepaired = inputPSDInfo.repaired;
info.inputPSDRelativeChange = norm(Q - Qoriginal / covarianceScale,'fro') / ...
    max(norm(Qoriginal / covarianceScale,'fro'),realmin);
info.psdEnforced = params.enforcePSD;
info.psdRepairApplied = psdInfo.applied || finalPSDInfo.repaired;
info.thetaScalePSD = psdInfo.thetaScale;
info.thetaRelativeChangePSD = psdInfo.thetaRelativeChange;
info.rangeProjectionRelativeChange = psdInfo.rangeProjectionRelativeChange;
info.mahalanobisSquaredBeforeShrink = psdInfo.mahalanobisSquaredBeforeShrink;
info.decompositionRelativeError = norm(Q - theta_outer - EeeT_est,'fro') / max(norm(Q,'fro'),realmin);
info.covarianceScale = covarianceScale;
info.requestedSpectralShift = requestedBeta * covarianceScale;
info.actualSpectralShift = beta * covarianceScale;
info.spectralShiftAdapted = beta ~= requestedBeta;
info.p2Underdetermined = p == 2;
info.physicalErrorIdentified = false;
info.estimator = 'box-L1 followed by PSD feasibility map';
if ~params.enforcePSD
    info.estimator = 'raw box-L1 diagnostic (PSD not enforced)';
end
info.params = params;

% Return physical covariance units; all numerical tolerances above were relative.
EeeT_est = EeeT_est * covarianceScale;
theta_est = theta_est * sqrt(covarianceScale);

if params.verbose
    fprintf('NC: best objective = %.6g, best start = %d/%d, min diag(E) = %.6g, min eig(E) = %.6g\n', ...
        info.objective, info.bestStart, info.nStarts, info.minErrorVariance, ...
        info.minErrorCovEigenvalue);
end

end


%% ============================================================
% Local function 0: Option reader
%% ============================================================
function value = getOpt(opts,name,defaultValue)

if isfield(opts,name) && ~isempty(opts.(name))
    value = opts.(name);
else
    value = defaultValue;
end

end


function validate_parameters(params)
nonnegative = {'betaFactor','tol','objTol','psdTol'};
positive = {'eta0','etaMin'};
integers = {'maxIters','minIters','patience','lineSearchMax','nStarts','nRandomStarts','randomSeed'};
for i = 1:numel(nonnegative)
    validateattributes(params.(nonnegative{i}), {'numeric'}, ...
        {'real','scalar','finite','nonnegative'}, 'NC', nonnegative{i});
end
for i = 1:numel(positive)
    validateattributes(params.(positive{i}), {'numeric'}, ...
        {'real','scalar','finite','positive'}, 'NC', positive{i});
end
for i = 1:numel(integers)
    validateattributes(params.(integers{i}), {'numeric'}, ...
        {'real','scalar','finite','integer','nonnegative'}, 'NC', integers{i});
end
if params.maxIters < 1 || params.lineSearchMax < 1 || params.nStarts < 1
    error('NC:InvalidOptions', 'maxIters, lineSearchMax, and nStarts must be positive.');
end
validateattributes(params.etaDecay, {'numeric'}, {'real','scalar','finite','>',0,'<',1}, 'NC','etaDecay');
validateattributes(params.enforcePSD, {'logical','numeric'}, {'scalar','binary'}, 'NC','enforcePSD');
validateattributes(params.repairInputPSD, {'logical','numeric'}, {'scalar','binary'}, 'NC','repairInputPSD');
end


%% ============================================================
% Local function 1: Initializations
%% ============================================================
function thetaStarts = build_theta_starts(Q,beta,sqrtDiagQ,params)

p = size(Q,1);

[vecs, vals] = eig(Q - beta * eye(p));
eigVals = diag(vals);
[eigVals, order] = sort(eigVals, 'descend');
vecs = vecs(:,order);

nSpectral = min(params.nStarts, p);
thetaStarts = zeros(p, nSpectral + params.nRandomStarts);

for i = 1:nSpectral
    theta = sqrt(max(eigVals(i),0)) * vecs(:,i);
    thetaStarts(:,i) = project_theta(theta,sqrtDiagQ);
end

if params.nRandomStarts > 0
    stream = RandStream('mt19937ar','Seed',params.randomSeed);
    for i = 1:params.nRandomStarts
        idx = nSpectral + i;
        theta = sqrtDiagQ .* (2*rand(stream,p,1) - 1);
        thetaStarts(:,idx) = project_theta(theta,sqrtDiagQ);
    end
end

% Remove duplicate zero starts except when every initialization is zero.
norms = vecnorm(thetaStarts);
if any(norms > 0)
    thetaStarts = thetaStarts(:, norms > 0);
end

end


%% ============================================================
% Local function 2: Projected L1 optimization
%% ============================================================
function [bestTheta,bestObj,info] = run_projected_l1(Q,theta0,sqrtDiagQ,offMask,params)

theta = project_theta(theta0,sqrtDiagQ);
obj_prev = offdiag_L1_objective(Q,theta,offMask);

bestTheta = theta;
bestObj = obj_prev;

eta = params.eta0 / max(norm(theta),1e-8);
noImprove = 0;
converged = false;

for k = 1:params.maxIters

    theta_outer = theta * theta';

    % Exact subgradient of sum_{i~=j}|Q_ij - theta_i*theta_j|.
    % Both (i,j) and (j,i) occur in the symmetric sum.
    S = sign(theta_outer - Q);
    S(~offMask) = 0;
    grad = 2 * S * theta;

    accepted = false;
    etaTry = eta;

    for j = 1:params.lineSearchMax
        theta_new = theta - etaTry * grad;
        theta_new = project_theta(theta_new,sqrtDiagQ);
        obj_new = offdiag_L1_objective(Q,theta_new,offMask);

        if obj_new <= obj_prev
            accepted = true;
            break;
        end

        etaTry = max(etaTry * params.etaDecay, params.etaMin);
    end

    if ~accepted
        obj_new = obj_prev;
        theta_new = theta;
    end

    eta = etaTry;

    if obj_new < bestObj - params.objTol * max(1,bestObj)
        bestObj = obj_new;
        bestTheta = theta_new;
        noImprove = 0;
    else
        noImprove = noImprove + 1;
    end

    if k >= params.minIters && noImprove >= params.patience
        converged = true;
        break;
    end

    theta = theta_new;
    obj_prev = obj_new;
end

info = struct();
info.initialObjective = offdiag_L1_objective(Q,theta0,offMask);
info.finalObjective = bestObj;
info.iterations = k;
info.converged = converged;

end


%% ============================================================
% Local function 3: Projection
%% ============================================================
function theta = project_theta(theta,sqrtDiagQ)

% Enforce |theta_i| <= sqrt(Q_ii)
theta = sign(theta) .* min(abs(theta),sqrtDiagQ);

end


%% ============================================================
% Local function 4: Off-diagonal L1 objective
%% ============================================================
function obj = offdiag_L1_objective(Q,theta,offMask)

R = theta*theta' - Q;
obj = sum(abs(R(offMask)));

end


%% ============================================================
% Local function 5: PSD handling
%% ============================================================
function [theta, info] = enforce_residual_psd(Q, theta, tol)

% For Q >= 0, Q - theta*theta' >= 0 iff theta is in range(Q) and
% theta' * pinv(Q) * theta <= 1. The numerical range uses a roundoff cutoff.
thetaBefore = theta;
info = struct('applied',false,'thetaScale',1, ...
    'rangeProjectionRelativeChange',0,'thetaRelativeChange',0, ...
    'mahalanobisSquaredBeforeShrink',0, ...
    'minEigenvalueBefore',min(eig(Q-theta*theta')), ...
    'minEigenvalueAfter',NaN);
[V, D] = eig(0.5 * (Q + Q'));
qEig = real(diag(D));
qScale = max(abs(qEig));
pos = qEig > size(Q,1) * eps * qScale;
if all(pos)
    thetaRot = V' * theta;
else
    thetaRot = V(:,pos)' * theta;
    theta = V(:,pos) * thetaRot;
end
info.rangeProjectionRelativeChange = norm(theta - thetaBefore) / max(norm(thetaBefore),realmin);
gamma = sum(thetaRot.^2 ./ qEig(pos));
info.mahalanobisSquaredBeforeShrink = gamma;
alpha = 1;
if gamma >= 1
    % A machine-level interior margin avoids an artificial negative eigenvalue.
    alpha = max(0, 1 - 16 * size(Q,1) * eps) / sqrt(gamma);
end
theta = alpha * theta;
info.thetaScale = alpha;
info.thetaRelativeChange = norm(theta - thetaBefore) / max(norm(thetaBefore),realmin);
info.applied = info.thetaRelativeChange > 16 * size(Q,1) * eps;
info.minEigenvalueAfter = min(eig(Q - theta*theta'));
if info.minEigenvalueAfter < -tol * max(norm(Q,'fro'),realmin)
    error('NC:PSDProjectionFailure','PSD feasibility map exceeded the relative numerical tolerance.');
end

end


function [A, info] = repair_psd_matrix(A, tol, allowRepair, labelText)

A = 0.5 * (A + A');
info = psd_matrix_info(A);

scaleA = max(norm(A, 'fro'),realmin);
absTol = tol * scaleA;
if info.minEigenvalueBefore >= 0
    info.minEigenvalueAfter = info.minEigenvalueBefore;
    return;
end

if ~allowRepair && info.minEigenvalueBefore < -absTol
    error('%s is not positive semidefinite. Minimum eigenvalue: %.6g.', ...
        labelText, info.minEigenvalueBefore);
elseif ~allowRepair
    info.minEigenvalueAfter = info.minEigenvalueBefore;
    return;
end

[V, D] = eig(A);
eigVals = real(diag(D));
eigVals(eigVals < 0) = 0;
A = V * diag(eigVals) * V';
A = 0.5 * (A + A');

info.repaired = true;
info.minEigenvalueAfter = min(eig(A));

end


function info = psd_matrix_info(A)

eigVals = eig(0.5 * (A + A'));
info = struct();
info.repaired = false;
info.minEigenvalueBefore = min(eigVals);
info.minEigenvalueAfter = info.minEigenvalueBefore;

end
