# MATLAB Implementation

This directory provides MATLAB scripts for reproducing SNR-based data merging experiments.

➡️ For theoretical background and method descriptions, see the top-level README.

---

## Files

### Main Scripts

- **mergingExample_the.m**  
  Theoretical example based on covariance matrices (no time series)

- **mergingExample_syn.m**  
  Synthetic time-series example with explicit signal and error generation

---

### Supporting Functions

- **WA.m**  
  Weighted averaging based on error covariance

- **SNRopt.m**  
  SNR-optimal merging weights

- **maxR.m**  
  Correlation-maximizing weights

- **SNRest.m**  
  Estimates scaling factor and noise-to-signal ratio

- **NC.m**  
  Estimates error covariance and correlation (generalized TC)

- **EeeTGEN.m**  
  Generates synthetic error covariance matrix

- **dataGEN.m**  
  Generates synthetic signal and error

---

## How to Run

### 1. Add path

```matlab
addpath(genpath(pwd))
