# MATLAB scripts

Run the following from the repository root.

```matlab
addpath('MATLAB scripts')
mergingExample_the
mergingExample_syn
TCvsNC
```

The two merging examples print MSE and squared correlation for WA, SNR-opt, maxR, and EW and display the comparisons. `TCvsNC.m` prints `sweep_summary` and plots the paired variance differences, median error-variance MSE, and failure rates across the error-correlation sweep.

- `NC.m` estimates the error covariance, latent loadings, squared correlations, and diagnostics. It replaces the historical `ECVest.m`.
- `SNRest.m` uses one NC fit and returns `[N, a, snr, fit]`. Existing two-output calls are supported. Supply a target signal variance as the second argument, or omit it for a standardized target.
- `SNRopt.m`, `WA.m`, and `maxR.m` compute merging weights. `TC.m` is the three-product baseline.
- `EeeTGEN.m` and `dataGEN.m` generate synthetic covariances and time series.

[Model, usage, numerical results, and references](../README.md)
