# SNR-opt (Python Implementation)

This repository provides a Python implementation accompanying the paper:

Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022)  
Rethinking Satellite Data Merging: From Averaging to SNR Optimization  
IEEE Transactions on Geoscience and Remote Sensing, 60, 4405215  
https://doi.org/10.1109/TGRS.2021.3107028

---

## Overview

This repository implements and compares multiple data merging methods for multi-source observations:

- Weighted Averaging (WA)
- SNR-optimal merging (SNRopt)
- Correlation-maximizing merging (maxR)
- Equal Weighting (EW)

It also provides tools for:
- Synthetic data generation
- Covariance-based parameter estimation
- Performance evaluation using MSE and R²

---

## Code Structure

The repository consists of one main notebook and one supporting module.

---

### Main Script

1. **SNRopt_modular.ipynb**  
   - Interactive example (time-series based)  
   - Demonstrates full workflow step-by-step  
   - Compares merging performance under true and estimated parameters  

---

### Supporting Module

1. **snropt_helpers.py**  
   Provides all core functions required for the workflow:

- **eeeT_gen(p, ecc, var_e=None, tol=1e-10)**  
  Generates synthetic error covariance matrix  

- **data_gen(n, p, ecc, snr_db, tol=1e-10)**  
  Generates synthetic signal and error such that:

$$
x = a y + e
$$

- **wa(EeeT)**  
  Computes merging weights using inverse error covariance  

- **snr_opt(N, a)**  
  Computes SNR-optimal weights by solving:

$$
(N + a a^T) u = a
$$

- **maxR(theta, Q)**  
  Computes weights that maximize correlation with the true signal  

- **snr_est(ExxT, Ey2, tol=1e-12)**  
  Estimates:
  - scaling factor $a$  
  - noise-to-signal ratio matrix $N$  

- **nc(covx, tol=1e-12)**  
  Estimates:
  - error covariance matrix  
  - data–truth correlation  

  Equivalent to Triple Collocation when $p = 3$,  
  but stable for any number of datasets  

- **evaluate_metrics(ExxT, Ey2, a, U)**  
  Evaluates MSE and R² for given merging weights  

---

## Signal + Noise Model

$$
x = a y + e
$$

$$
Q = E[xx^T] = a a^T E[y^2] + E[ee^T]
$$

---

## SNR-based Decomposition

$$
C = \frac{Q}{E[y^2]}
$$

$$
C = a a^T + N
$$

---

## SNR-optimal Weights

$$
u = (N + a a^T)^{-1} a
$$

---

## Typical Workflow

1. Generate synthetic data  
2. Compute covariance matrices  
3. Estimate parameters (NC, SNRest)  
4. Compute merging weights (WA, SNRopt, maxR, EW)  
5. Evaluate performance using MSE and R²  

---

## Requirements

Install required packages:

```bash
pip install numpy matplotlib jupyter
