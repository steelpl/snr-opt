# Jupyter Notebook

Open a terminal in this directory and run the following. Open the notebook and select **Run All**.

```bash
python -m pip install numpy matplotlib notebook
python -m notebook SNRopt_modular.ipynb
```

[SNRopt_modular.ipynb](SNRopt_modular.ipynb) generates collocated synthetic time series, compares four merging methods using true and NC-estimated statistics, and displays MSE and squared correlation. It imports the estimator and weight functions from [snropt_helpers.py](snropt_helpers.py).

- `nc(Q)` returns error covariance, latent loadings, and squared correlations. Pass `return_info=True` for NC diagnostics.
- `snr_est(Q, Ey2)` returns normalized noise covariance and target-scale loadings. Pass `return_fit=True` to also obtain the per-product SNR and the fit dictionary. Omit `Ey2` for a unit-variance target.
- `snr_opt(N, a)`, `wa(E)`, and `maxR(theta, Q)` compute merging weights.

[Model, usage, numerical results, and references](../README.md)
