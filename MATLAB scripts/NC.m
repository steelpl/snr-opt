function [EeeT_est, theta_est, rho2_est, info] = NC(covx, opts)
% ============================================================
% N-Tuple Collocation (NC)
% Adapted from NC-public (be461e3, MIT) for this study; see README.md.
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
% Revised estimator for the three-merging-method comparison. Returned estimates
% use an admissible exact solution for three products. Other fits use
% dimensionless off-diagonal L1 loss; feasible residual candidates are
% compared after PSD projection. For four or more products this is a
% multistart numerical estimate, without a global optimality certificate.
% info.objective retains the original covariance-unit residual for callers.
% info.standardizedObjective and info.relativeResidual describe the fitted
% dimensionless loss. FitClass is exact/bound/sign for nonzero triples,
% degenerate for triples with zero covariances, multi for larger inputs,
% and underdetermined for two products. Polarities are anchored to the
% first non-negligible standardized loading.
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

%% Work in dimensionless product units.
Qoriginal = 0.5*double(covx) + 0.5*double(covx');
covarianceScale=max(abs(Qoriginal(:)));
if covarianceScale==0, covarianceScale=1; end
Qbase=Qoriginal;
inputRepaired=false;
% An invalid covariance can contain a non-positive variance with nonzero
% cross terms. Retain the original opt-in repair behavior for these inputs.
zeroVar=diag(Qbase)<=0;
if any(diag(Qbase)<0) || any(any(Qbase(zeroVar,:)~=0))
    [Qbase,initial]=repair_psd_matrix(Qbase,params.psdTol,params.repairInputPSD,'Input covariance matrix');
    inputRepaired=initial.repaired;
end
scale=sqrt(max(0,diag(Qbase)));
divisor=scale; divisor(divisor==0)=1;
R=(Qbase./divisor)./divisor';
[R,check]=repair_psd_matrix(R,params.psdTol,params.repairInputPSD,'Input product correlation matrix');
inputRepaired=inputRepaired||check.repaired;
Q=(R.*divisor).*divisor';
diagR=max(0,diag(R));
bounds=sqrt(diagR);
offMask=~eye(p);
beta=min(params.betaFactor*min(diagR),0.95*max(eig(R)));
fitClass='multi'; exact=false;
tri=[];
if p==3
    qr=[R(1,2);R(1,3);R(2,3)];
    if any(qr==0)
        fitClass='degenerate';
    elseif prod(sign(qr))<0
        fitClass='sign';
    else
        sq=[qr(1)*(qr(2)/qr(3)); qr(1)*(qr(3)/qr(2)); qr(2)*(qr(3)/qr(1))];
        tri=sqrt(sq).*[1;sign(qr(1));sign(qr(2))];
        exact=all(sq<=diagR+16*eps(max(1,max(diagR))));
        if exact, fitClass='exact'; else, fitClass='bound'; end
    end
elseif p==2
    fitClass='underdetermined';
end
% Keep the objective and the PSD map in the same dimensionless units. The
% rounded correlation is only used to generate search candidates for p>=4.
% Selection and all reported diagnostics use the original correlation.
if exact
    u=min(max(tri,-bounds),bounds);
    raw=u;
    bestStart=1;
    startInfo=struct('initialObjective',0,'finalObjective',offdiag_L1_objective(R,u,offMask),'iterations',0,'converged',true);
    psdInfo=no_repair();
else
    searchR=R;
    if p>=4
        searchR=round(R,13);
    end
    anchor=R./max(bounds',realmin);
    if p==3 && ~any(diagR==0) && ~strcmp(fitClass,'degenerate')
        % For a three-product, nonzero, inadmissible system, a unit-bound
        % anchor can attain the global box-L1 minimum. Each such candidate
        % has a PSD residual (the corresponding Schur complement). A
        % bound violation leaves one underestimated edge; a sign conflict
        % must mismatch at least one edge. Changing either matched edge by
        % |d| changes their product by at most |d|. Evaluate all three.
        candidates=anchor;
    else
        spectral=build_theta_starts(searchR,beta,bounds,params);
        candidates=[anchor,spectral];
        for j=1:size(spectral,2)
            legacy=run_projected_l1(searchR,spectral(:,j),bounds,offMask,params);
            candidates(:,end+1)=legacy; %#ok<AGROW>
        end
    end
    if ~isempty(tri), candidates(:,end+1)=project_theta(tri,bounds); end
    count=size(candidates,2);
    startInfo=repmat(struct('initialObjective',0,'finalObjective',0,'iterations',0,'converged',true),count,1);
    best=inf; bestStart=1;
    for j=1:count
        candidate=candidates(:,j);
        initialObj=offdiag_L1_objective(R,candidate,offMask);
        if p>=4 || strcmp(fitClass,'degenerate') || p==2
            [trial,it]=coordinate_medians(searchR,candidate,bounds,params);
            proposals=[candidate,trial];
        else
            it=0;
            proposals=candidate;
        end
        finalObj=Inf;
        for k=1:size(proposals,2)
            rawTrial=proposals(:,k); uTrial=rawTrial; map=no_repair();
            if params.enforcePSD && min(eig(R-uTrial*uTrial')) < -16*p*eps(max(1,norm(R,'fro')))
                [uTrial,map]=enforce_residual_psd(R,uTrial,params.psdTol);
            end
            obj=offdiag_L1_objective(R,uTrial,offMask);
            finalObj=min(finalObj,obj);
            if obj < best-32*p*eps(max(1,obj)) || isinf(best)
                best=obj; u=uTrial; raw=rawTrial; psdInfo=map; bestStart=j;
            end
        end
        startInfo(j)=struct('initialObjective',initialObj,'finalObjective',finalObj,'iterations',it,'converged',it<params.maxIters);
    end
end
% The overall sign is arbitrary. Anchor to the first non-negligible
% standardized loading so multiplying a product by 100 cannot flip the score.
anchor=find(abs(u)>32*eps*max(abs(u)),1);
if ~isempty(anchor) && u(anchor)<0, u=-u; raw=-raw; end
theta_est=divisor.*u;
EeeT_est=Q-theta_est*theta_est';
EeeT_est=(EeeT_est+EeeT_est')/2;
rho2_est=nan(p,1); positive=diagR>0;
rho2_est(positive)=min(1,max(0,u(positive).^2./diagR(positive)));
offTotal=sum(abs(R(offMask)));
rawTotal=sum(abs(Q(offMask)));
rawObj=offdiag_L1_objective(Q,theta_est,offMask);
rawBefore=offdiag_L1_objective(Q,divisor.*raw,offMask);
stdObj=offdiag_L1_objective(R,u,offMask);
residual=R-u*u';
info=struct();
info.FitClass=fitClass;
info.standardizedObjective=stdObj;
info.standardizedObjectiveBeforePSD=offdiag_L1_objective(R,raw,offMask);
info.relativeResidual=stdObj/max(realmin,offTotal);
info.rawRelativeResidual=rawObj/max(realmin,rawTotal);
info.boundActive=positive & abs(abs(raw)-bounds)<=1e-8.*bounds;
info.boundActiveAfterPSD=positive & abs(abs(u)-bounds)<=1e-8.*bounds;
info.negativeLoading=u < -1e-10*max(abs(u));
weight=pinv(R)*u;
info.negativeWeight=weight < -1e-10*max(abs(weight));
info.objective=rawObj;
info.objectiveBeforePSD=rawBefore;
info.objectiveChangePSD=rawObj-rawBefore;
info.objectiveChangePSDRelativeToQ=(rawObj-rawBefore)/max(norm(Q,'fro'),realmin);
info.bestStart=bestStart;
info.nStarts=numel(startInfo);
info.startInfo=startInfo; % start objectives are now in correlation units.
info.minErrorVariance=min(diag(EeeT_est));
info.minErrorCovEigenvalue=min(eig(EeeT_est));
info.maxRho2=max(rho2_est);
info.feasibleDiagonal=all(diag(residual)>=-params.psdTol);
info.feasiblePSD=min(eig(residual))>=-params.psdTol*max(norm(R,'fro'),realmin);
info.inputMinEigenvalue=min(eig(Qoriginal));
info.inputPSDRepaired=inputRepaired;
info.inputPSDRelativeChange=norm(Q-Qoriginal,'fro')/max(norm(Qoriginal,'fro'),realmin);
info.psdEnforced=params.enforcePSD;
info.psdRepairApplied=psdInfo.applied;
info.thetaScalePSD=psdInfo.thetaScale;
info.thetaRelativeChangePSD=psdInfo.thetaRelativeChange;
info.rangeProjectionRelativeChange=psdInfo.rangeProjectionRelativeChange;
info.mahalanobisSquaredBeforeShrink=psdInfo.mahalanobisSquaredBeforeShrink;
info.decompositionRelativeError=norm(Q-theta_est*theta_est'-EeeT_est,'fro')/max(norm(Q,'fro'),realmin);
info.covarianceScale=covarianceScale;
info.requestedSpectralShift=params.betaFactor*min(diagR);
info.actualSpectralShift=beta;
info.spectralShiftAdapted=beta~=info.requestedSpectralShift;
info.p2Underdetermined=p==2;
info.physicalErrorIdentified=false;
info.estimator='Exact three-product fit or dimensionless off-diagonal L1 fit';
info.params=params;
if params.verbose
    fprintf('NC: %s, correlation-unit relative residual %.5g, minimum error eigenvalue %.5g\n',fitClass,info.relativeResidual,info.minErrorCovEigenvalue);
end
end

function result=no_repair()
result=struct('applied',false,'thetaScale',1,'rangeProjectionRelativeChange',0,'thetaRelativeChange',0,'mahalanobisSquaredBeforeShrink',NaN);
end

function [u,iter]=coordinate_medians(R,u,bounds,params)
% A coordinatewise minimum is the weighted median of R_ij/u_j with weights
% |u_j|. Choose the projection of the current value into a tied interval.
p=length(u);
for iter=1:params.maxIters
    last=u;
    for i=1:p
        keep=u~=0; keep(i)=false;
        if ~any(keep), continue; end
        [values,order]=sort(R(keep,i)./u(keep));
        weights=abs(u(keep)); cumulative=cumsum(weights(order)); half=cumulative(end)/2;
        tol=16*eps(cumulative(end));
        left=find(cumulative>=half-tol,1);
        right=find(cumulative>=half+tol,1);
        if isempty(right),right=numel(values);end
        u(i)=max(-bounds(i),min(bounds(i),max(values(left),min(values(right),u(i)))));
    end
    if norm(u-last,Inf)<=max(32*eps,params.tol), break; end
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
