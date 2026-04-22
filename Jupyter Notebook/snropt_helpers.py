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

    EeeT = np.diag(var_e)
    for i in range(p):
        for j in range(i + 1, p):
            val = np.sqrt(var_e[i] * var_e[j]) * ecc
            EeeT[i, j] = val
            EeeT[j, i] = val

    EeeT = 0.5 * (EeeT + EeeT.T)
    vals, vecs = np.linalg.eigh(EeeT)
    vals = np.maximum(vals, tol)
    return vecs @ np.diag(vals) @ vecs.T


def data_gen(n: int, p: int, ecc: float, snr_db: float, tol: float = 1e-10):
    if n <= 1:
        raise ValueError("n must be greater than 1.")
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

    ye = ye @ np.linalg.inv(Rc)
    ye = ye @ Rm

    y = ye[:, 0]
    e = ye[:, 1:]
    return y, e


def wa(EeeT: np.ndarray) -> np.ndarray:
    EeeT = 0.5 * (EeeT + EeeT.T)
    eta = np.ones((EeeT.shape[0], 1))
    tmp = np.linalg.solve(EeeT, eta)
    u = tmp / float(eta.T @ tmp)
    return u.ravel()


def maxR(theta: np.ndarray, Q: np.ndarray) -> np.ndarray:
    theta = np.asarray(theta, dtype=float).reshape(-1, 1)
    Q = 0.5 * (Q + Q.T)
    u = np.linalg.solve(Q, theta).ravel()
    s = u.sum()
    if abs(s) > 1e-12:
        u = u / s
    return u


def snr_opt(N: np.ndarray, a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float).reshape(-1, 1)
    N = 0.5 * (N + N.T)
    u = np.linalg.solve(N + a @ a.T, a)
    return u.ravel()


def snr_est(ExxT: np.ndarray, Ey2: float, tol: float = 1e-12):
    ExxT = np.asarray(ExxT, dtype=float)

    if ExxT.ndim != 2 or ExxT.shape[0] != ExxT.shape[1]:
        raise ValueError("ExxT must be a square matrix.")
    if Ey2 <= 0:
        raise ValueError("Ey2 must be positive.")

    Q = 0.5 * (ExxT + ExxT.T)
    C = Q / Ey2
    p = C.shape[0]

    diagC = np.diag(C)
    if np.any(diagC <= 0):
        raise ValueError("Normalized covariance has non-positive diagonal entries.")

    sqrtDiagC = np.sqrt(diagC)

    beta = 0.5 * np.min(diagC)
    step = 0.01
    iters = 2000

    vals, vecs = np.linalg.eigh(C - beta * np.eye(p))
    idx = np.argmax(vals)
    a_est = np.sqrt(max(vals[idx], 0.0)) * vecs[:, idx]

    # same role as MATLAB sign-consistency step
    a_est = np.abs(a_est)

    step = step / max(np.linalg.norm(a_est), 1e-8)

    mask = np.ones((p, p)) - np.eye(p)

    for _ in range(iters):
        A = np.outer(a_est, a_est)
        S = np.sign(A - C)
        S[np.diag_indices(p)] = 0.0
        grad = S @ a_est
        a_est = a_est - step * grad
        a_est = np.sign(a_est) * np.minimum(np.abs(a_est), sqrtDiagC)

    N_est = C - np.outer(a_est, a_est)
    N_est = 0.5 * (N_est + N_est.T)

    d = np.diag(N_est).copy()
    d[(d < 0) & (d > -tol)] = 0.0
    N_est[np.diag_indices(p)] = d

    return N_est, a_est


def nc(covx: np.ndarray, tol: float = 1e-12):
    Q = 0.5 * (covx + covx.T)
    p = Q.shape[0]
    diagQ = np.diag(Q)
    if np.any(diagQ < 0):
        raise ValueError("Input covariance has negative diagonal entries.")
    sqrtDiagQ = np.sqrt(np.maximum(diagQ, 0.0))

    beta = 0.5 * np.min(diagQ)
    step = 0.01
    iters = 2000
    mask = np.ones((p, p)) - np.eye(p)

    vals, vecs = np.linalg.eigh(Q - beta * np.eye(p))
    idx = np.argmax(vals)
    theta_est = np.sqrt(max(vals[idx], 0.0)) * vecs[:, idx]
    theta_est = np.sign(theta_est) * np.abs(theta_est)
    step = step / max(np.linalg.norm(theta_est), 1e-8)

    for _ in range(iters):
        theta_outer = np.outer(theta_est, theta_est)
        grad = (mask * np.sign(theta_outer - Q)) @ theta_est
        theta_est = theta_est - step * grad
        theta_est = np.sign(theta_est) * np.minimum(np.abs(theta_est), sqrtDiagQ)

    EeeT_est = Q - np.outer(theta_est, theta_est)
    EeeT_est = 0.5 * (EeeT_est + EeeT_est.T)
    d = np.diag(EeeT_est).copy()
    d[(d < 0) & (d > -tol)] = 0.0
    EeeT_est[np.diag_indices(p)] = d

    rho2_est = (theta_est ** 2) / np.maximum(diagQ, tol)
    rho2_est[(rho2_est < 0) & (rho2_est > -tol)] = 0.0
    rho2_est[(rho2_est > 1) & (rho2_est < 1 + tol)] = 1.0
    return EeeT_est, theta_est, rho2_est


def evaluate_metrics(ExxT: np.ndarray, Ey2: float, a: np.ndarray, U: np.ndarray):
    a = np.asarray(a, dtype=float).reshape(-1, 1)
    U = np.asarray(U, dtype=float)
    if U.ndim == 1:
        U = U.reshape(-1, 1)

    mse = np.diag(U.T @ ExxT @ U - 2 * Ey2 * (U.T @ a) + Ey2)
    num = np.diag(U.T @ (a @ a.T) @ U)
    den = np.maximum(np.diag(U.T @ ExxT @ U), 1e-12)
    r2 = Ey2 * num / den
    return mse, r2
