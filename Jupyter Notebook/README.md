# SNRopt Modular Example

This repository contains a modular Python implementation of synthetic data generation, covariance-based parameter estimation, and data-merging methods for comparing WA, SNRopt, and maxR.

## Files

- `SNRopt_modular.ipynb`  
  Main practice notebook for running the workflow step by step.

- `snropt_helpers.py`  
  Helper module containing the core functions used by the notebook.

## Included Functions

The helper module provides the following functions:

- `eeeT_gen(p, ecc, var_e=None, tol=1e-10)`  
  Generates a synthetic error covariance matrix.

- `data_gen(n, p, ecc, snr_db, tol=1e-10)`  
  Generates synthetic zero-mean signal and error time series.

- `wa(EeeT)`  
  Computes weighted-average merging weights from the error covariance matrix.

- `maxR(theta, Q)`  
  Computes correlation-maximizing weights.

- `snr_opt(N, a)`  
  Computes SNR-optimal merging weights from the normalized signal-plus-noise representation.

- `snr_est(ExxT, Ey2, tol=1e-12)`  
  Estimates the normalized scaling vector and noise-to-signal matrix.

- `nc(covx, tol=1e-12)`  
  Estimates error covariance, signal scaling, and squared correlation using N-Tuple Collocation.

- `evaluate_metrics(ExxT, Ey2, a, U)`  
  Evaluates MSE and R² for one or more weight vectors.

## Requirements

Install the required Python packages before running the notebook:

```bash
pip install numpy matplotlib jupyter
