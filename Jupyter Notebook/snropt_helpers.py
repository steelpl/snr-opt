"""Python helpers for SNR-opt and N-tuple collocation.

SNR-opt: Kim et al. (2022), https://doi.org/10.1109/TGRS.2021.3107028
NC: Pham & Kim (2027), https://doi.org/10.1016/j.jhydrol.2026.136749
NC algorithm adapted from steelpl/NC-public at be461e3 (MIT).
"""

import numpy as np


def eeeT_gen(p: int, ecc: float, var_e=None, tol: float = 1e-10) -> np.ndarray:
    if p < 1:
        raise ValueError("p must be a positive integer.")
    ecc = float(np.clip(ecc, 0.0, 1.0))
    if var_e is None:
        var_e = np.random.rand(p)
    var_e = np.asarray(var_e, dtype=float).reshape(-1)
    if len(var_e) != p:
        raise ValueError("var_e must have length p.")
    if np.any(var_e < 0):
        raise ValueError("var_e must be non-negative.")

    sd = np.sqrt(var_e)
    EeeT = (1 - ecc) * np.diag(var_e) + ecc * np.outer(sd, sd)
    vals, vecs = np.linalg.eigh(EeeT)
    vals = np.maximum(vals, tol)
    return vecs @ np.diag(vals) @ vecs.T


def data_gen(n: int, p: int, ecc: float, snr_db: float, tol: float = 1e-10):
    if n <= 1:
        raise ValueError("n must be greater than 1.")
    if not np.isfinite(snr_db):
        raise ValueError("snr_db must be finite.")
    snr = 10 ** (snr_db / 10.0)

    EeeT = eeeT_gen(p, ecc, tol=tol)
    Ey2 = np.mean(np.diag(EeeT) * snr)

    m = np.zeros((p + 1, p + 1), dtype=float)
    m[0, 0] = Ey2
    m[1:, 1:] = EeeT
    m = 0.5 * (m + m.T)

    vals, vecs = np.linalg.eigh(m)
    vals = np.maximum(vals, tol)
    m = vecs @ np.diag(vals) @ vecs.T

    ye = np.random.randn(n, p + 1)
    ye = ye - ye.mean(axis=0, keepdims=True)

    C = np.cov(ye, rowvar=False)
    C = 0.5 * (C + C.T)
    vals_c, vecs_c = np.linalg.eigh(C)
    vals_c = np.maximum(vals_c, tol)
    C = vecs_c @ np.diag(vals_c) @ vecs_c.T

    # MATLAB chol returns upper triangular R such that A = R.T @ R
    Rc = np.linalg.cholesky(C).T
    Rm = np.linalg.cholesky(m).T

    ye = np.linalg.solve(Rc.T, ye.T).T @ Rm

    y = ye[:, 0]
    e = ye[:, 1:]
    return y, e


def _symmetric(matrix, min_size=1):
    matrix = np.asarray(matrix)
    if (
        matrix.ndim != 2
        or matrix.shape[0] != matrix.shape[1]
        or matrix.shape[0] < min_size
        or np.iscomplexobj(matrix)
    ):
        raise ValueError(f"Expected a real square matrix of size >= {min_size}.")
    matrix = matrix.astype(float, copy=False)
    if not np.all(np.isfinite(matrix)):
        raise ValueError("Expected a finite matrix.")
    return 0.5 * matrix + 0.5 * matrix.T


def wa(EeeT: np.ndarray, covariance_scale=0.0) -> np.ndarray:
    """Minimize error variance subject to sum(u) = 1 for common-scale products.

    Use the observation covariance scale for the roundoff threshold when
    EeeT is computed by subtraction in NC.
    """
    values, vectors = np.linalg.eigh(_symmetric(EeeT))
    p = values.size
    eps = np.finfo(float).eps
    if not np.isscalar(covariance_scale) or not np.isfinite(covariance_scale) or covariance_scale < 0:
        raise ValueError("covariance_scale must be a finite nonnegative scalar.")
    tolerance = 64 * p * eps * max(np.max(np.abs(values)), covariance_scale)
    if np.any(values < -tolerance):
        raise ValueError("Error covariance must be positive semidefinite.")
    ones = np.ones(p)
    positive = values > tolerance
    null_vectors = vectors[:, ~positive]
    null_component = null_vectors.T @ ones
    if np.linalg.norm(null_component) > 64 * p * eps * np.linalg.norm(ones):
        direction = null_vectors @ null_component
    else:
        basis = vectors[:, positive]
        eigenvalues = values[positive]
        direction = basis @ ((basis.T @ ones) * (eigenvalues.min() / eigenvalues))
    return direction / direction.sum()


def maxR(theta: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Correlation-maximizing weights, normalized to sum to one when stable."""
    Q = _symmetric(Q)
    theta = np.asarray(theta, dtype=float).reshape(-1)
    u = np.linalg.pinv(Q, rcond=len(Q) * np.finfo(float).eps) @ theta
    total = u.sum()
    if abs(total) > 1e-12 * np.linalg.norm(u, ord=1):
        u = u / total
    return u


def snr_opt(N: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Linear minimum-MSE weights for the normalized signal/noise model."""
    a = np.asarray(a, dtype=float).reshape(-1)
    C = _symmetric(N) + np.outer(a, a)
    return np.linalg.pinv(C, rcond=len(C) * np.finfo(float).eps) @ a


def _projected_l1(Q, theta, bound, off_diagonal, params):
    def objective(candidate):
        residual = np.outer(candidate, candidate) - Q
        return np.abs(residual[off_diagonal]).sum()

    previous = best_loss = initial_loss = objective(theta)
    best = theta
    step = params["eta0"] / max(np.linalg.norm(theta), 1e-8)
    stagnant = 0
    converged = False
    for iteration in range(1, params["maxIters"] + 1):
        signs = np.sign(np.outer(theta, theta) - Q)
        signs[~off_diagonal] = 0
        gradient = 2 * (signs @ theta)
        for _ in range(params["lineSearchMax"]):
            candidate = np.clip(theta - step * gradient, -bound, bound)
            loss = objective(candidate)
            if loss <= previous:
                break
            step = max(step * params["etaDecay"], params["etaMin"])
        else:
            candidate, loss = theta, previous

        if loss < best_loss - params["objTol"] * max(1, best_loss):
            best, best_loss, stagnant = candidate, loss, 0
        else:
            stagnant += 1
        if iteration >= params["minIters"] and stagnant >= params["patience"]:
            converged = True
            break
        theta, previous = candidate, loss

    return best, best_loss, {
        "initialObjective": initial_loss,
        "finalObjective": best_loss,
        "iterations": iteration,
        "converged": converged,
    }


def nc(covx: np.ndarray, tol: float = 1e-12, *, opts=None, return_info=False):
    """Fit box-constrained off-diagonal L1, then enforce residual PSD by default.

    Returns (error_cov, theta, rho2) or (error_cov, theta, rho2, info).
    The options use the same names and defaults as MATLAB NC.m.
    """
    Q = _symmetric(covx, min_size=2)
    p = len(Q)
    params = {
        "betaFactor": 0.5, "eta0": 0.005, "maxIters": 2000, "minIters": 20,
        "patience": 30, "tol": tol, "objTol": 1e-10, "etaDecay": 0.7,
        "etaMin": 1e-10, "lineSearchMax": 12, "nStarts": min(4, p),
        "nRandomStarts": 0, "randomSeed": 1, "enforcePSD": True,
        "repairInputPSD": True, "psdTol": 1e-10, "verbose": False,
    }
    if opts is not None:
        unknown = opts.keys() - params.keys()
        if unknown:
            raise ValueError(f"Unknown NC options: {sorted(unknown)}")
        params.update(opts)
    for name in ("enforcePSD", "repairInputPSD", "verbose"):
        if not np.isscalar(params[name]) or params[name] not in (0, 1):
            raise ValueError(f"{name} must be boolean.")
    for name in ("maxIters", "minIters", "patience", "lineSearchMax",
                 "nStarts", "nRandomStarts", "randomSeed"):
        value = params[name]
        if not np.isfinite(value) or value < 0 or value != int(value):
            raise ValueError(f"{name} must be a nonnegative integer.")
        params[name] = int(value)
    if not params["maxIters"] or not params["lineSearchMax"]:
        raise ValueError("maxIters and lineSearchMax must be positive.")
    for name in ("betaFactor", "tol", "objTol", "psdTol", "eta0", "etaMin", "etaDecay"):
        if not np.isfinite(params[name]) or params[name] < 0:
            raise ValueError(f"{name} must be finite and nonnegative.")
    if not (0 < params["etaDecay"] < 1) or not params["eta0"] or not params["etaMin"]:
        raise ValueError("Require positive steps and 0 < etaDecay < 1.")

    scale = np.max(np.abs(Q)) or 1.0
    Q = Q / scale
    original = Q
    values, vectors = np.linalg.eigh(Q)
    input_min = values[0]
    input_repaired = input_min < 0 and params["repairInputPSD"]
    if input_repaired:
        Q = (vectors * np.maximum(values, 0)) @ vectors.T
        Q = (Q + Q.T) / 2
    elif input_min < -params["psdTol"] * np.linalg.norm(Q):
        raise ValueError("Input covariance is not positive semidefinite.")
    variance = np.diag(Q).copy()
    if np.any(variance < -params["tol"]):
        raise ValueError("Input covariance has negative diagonal entries.")
    variance = np.maximum(variance, 0)
    Q = Q.copy()
    Q[np.diag_indices(p)] = variance
    bound = np.sqrt(variance)

    requested_shift = params["betaFactor"] * variance.min()
    shift = min(max(requested_shift, 0), 0.95 * max(np.linalg.eigvalsh(Q)[-1], 0))
    values, vectors = np.linalg.eigh(Q - shift * np.eye(p))
    n_starts = min(p, max(1, params["nStarts"]))
    starts = [
        np.clip(np.sqrt(max(values[i], 0)) * vectors[:, i], -bound, bound)
        for i in range(p - 1, p - n_starts - 1, -1)
    ]
    if params["nRandomStarts"]:
        rng = np.random.RandomState(params["randomSeed"])
        starts.extend(rng.uniform(-1, 1, (params["nRandomStarts"], p)) * bound)
    starts = [start for start in starts if np.any(start)] or starts
    off_diagonal = ~np.eye(p, dtype=bool)
    fits = [_projected_l1(Q, start, bound, off_diagonal, params) for start in starts]
    best = min(range(len(fits)), key=lambda i: fits[i][1])
    theta, loss_before, _ = fits[best]

    theta_before = theta
    shrink = 1.0
    gamma = np.nan
    range_change = 0.0
    eps = np.finfo(float).eps
    tiny = np.finfo(float).tiny
    if params["enforcePSD"]:
        values, vectors = np.linalg.eigh(Q)
        positive = values > p * eps * np.max(np.abs(values))
        basis = vectors[:, positive]
        rotated = basis.T @ theta
        if not positive.all():
            theta = basis @ rotated
        range_change = np.linalg.norm(theta - theta_before) / max(np.linalg.norm(theta_before), tiny)
        gamma = np.sum(rotated**2 / values[positive])
        if gamma >= 1:
            shrink = max(0, 1 - 16 * p * eps) / np.sqrt(gamma)
            theta = shrink * theta

    signal = np.outer(theta, theta)
    EeeT = Q - signal
    EeeT = (EeeT + EeeT.T) / 2
    min_error_eigenvalue = np.linalg.eigvalsh(EeeT)[0]
    feasible = min_error_eigenvalue >= -params["psdTol"] * max(np.linalg.norm(Q), tiny)
    if params["enforcePSD"] and not feasible:
        raise ValueError("NC residual exceeded the PSD numerical tolerance.")
    rho2 = np.full(p, np.nan)
    np.divide(theta**2, variance, out=rho2, where=variance > 0)
    rho2[(rho2 > 1) & (rho2 < 1 + params["tol"])] = 1.0

    relative_change = np.linalg.norm(theta - theta_before) / max(np.linalg.norm(theta_before), tiny)
    loss_after = np.abs(EeeT[off_diagonal]).sum()
    if params["verbose"]:
        print(f"NC: objective={loss_after * scale:.6g}, best start={best + 1}/{len(fits)}, "
              f"min eig(E)={min_error_eigenvalue * scale:.6g}")
    result = EeeT * scale, theta * np.sqrt(scale), rho2
    if not return_info:
        return result

    start_info = [
        dict(run, initialObjective=run["initialObjective"] * scale,
             finalObjective=run["finalObjective"] * scale)
        for _, _, run in fits
    ]
    info = {
        "objective": loss_after * scale,
        "objectiveBeforePSD": loss_before * scale,
        "objectiveChangePSD": (loss_after - loss_before) * scale,
        "bestStart": best + 1,
        "nStarts": len(fits),
        "startInfo": start_info,
        "minErrorVariance": np.diag(EeeT).min() * scale,
        "minErrorCovEigenvalue": min_error_eigenvalue * scale,
        "feasiblePSD": bool(feasible),
        "inputMinEigenvalue": input_min * scale,
        "inputPSDRepaired": bool(input_repaired),
        "inputPSDRelativeChange": np.linalg.norm(Q - original) / max(np.linalg.norm(original), tiny),
        "psdEnforced": bool(params["enforcePSD"]),
        "psdRepairApplied": bool(relative_change > 16 * p * eps),
        "thetaScalePSD": shrink,
        "thetaRelativeChangePSD": relative_change,
        "rangeProjectionRelativeChange": range_change,
        "mahalanobisSquaredBeforeShrink": gamma,
        "covarianceScale": scale,
        "p2Underdetermined": p == 2,
        "physicalErrorIdentified": False,
        "params": params,
    }
    return (*result, info)


def snr_est(ExxT: np.ndarray, Ey2: float = 1.0, tol: float = 1e-12,
            *, opts=None, return_fit=False):
    """Return (N, a) or (N, a, snr, fit) from one NC fit."""
    if not np.isscalar(Ey2) or np.iscomplexobj(Ey2) or not np.isfinite(Ey2) or Ey2 <= 0:
        raise ValueError("Ey2 must be a finite, positive scalar.")
    error_cov, theta, rho2, info = nc(ExxT, tol, opts=opts, return_info=True)
    if theta[np.argmax(np.abs(theta))] < 0:
        theta = -theta
    N, a = error_cov / Ey2, theta / np.sqrt(Ey2)
    if not return_fit:
        return N, a
    signal, noise = theta**2, np.maximum(np.diag(error_cov), 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        snr = signal / noise
    fit = dict(errorCov=error_cov, theta=theta, rho2=rho2, info=info,
               covariance=error_cov + np.outer(theta, theta))
    return N, a, snr, fit


def evaluate_metrics(ExxT: np.ndarray, Ey2: float, a: np.ndarray, U: np.ndarray):
    """Evaluate covariance-based MSE and squared correlation per weight column."""
    a = np.asarray(a, dtype=float).reshape(-1)
    U = np.asarray(U, dtype=float)
    if U.ndim == 1:
        U = U[:, None]
    merged_variance = np.sum(U * (ExxT @ U), axis=0)
    merged_gain = U.T @ a
    mse = merged_variance - 2 * Ey2 * merged_gain + Ey2
    r2 = np.full_like(merged_variance, np.nan)
    np.divide(Ey2 * merged_gain**2, merged_variance, out=r2, where=merged_variance > 0)
    return mse, r2
