# SNR-opt (MATLAB Implementation)

This repository provides MATLAB codes accompanying the paper:

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

The repository consists of two main scripts and seven supporting functions.

---

### Main Scripts

1. **mergingExample_the.m**  
   - Theoretical example (no time series)  
   - Uses analytical covariance matrices  
   - Demonstrates ideal performance under known parameters  

2. **mergingExample_syn.m**  
   - Synthetic time-series example  
   - Generates signal and error explicitly  
   - Demonstrates practical performance under estimated parameters  

---

### Supporting Functions

1. **WA.m**  
   Computes merging weights using inverse error covariance  

2. **SNRopt.m**  
   Computes SNR-optimal weights by solving:

$$
(N + a a^T) u = a
$$

3. **maxR.m**  
   Computes weights that maximize correlation with the true signal  

4. **SNRest.m**  
   Estimates:
   - scaling factor $a$  
   - noise-to-signal ratio matrix $N$  

5. **NC.m**  
   Estimates:
   - error covariance matrix  
   - data–truth correlation  

   Equivalent to Triple Collocation when $p = 3$,  
   but stable for any number of datasets  

6. **EeeTGEN.m**  
   Generates synthetic error covariance matrix  

7. **dataGEN.m**  
   Generates synthetic signal and error such that:

$$
x = a y + e
$$

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

1. Generate synthetic data or covariance matrix  
2. Compute merging weights (WA, SNRopt, maxR, EW)  
3. Evaluate performance using MSE and R²  
4. Compare true vs estimated parameter cases  

---

## Important Notes

- SNRopt is optimal when parameters are known  
- Performance degrades when parameters are estimated  
- WA is more robust but suboptimal  
- NC provides a stable alternative to TC  

---

## Purpose

This repository is intended to:

- Reproduce results from the paper  
- Provide transparent implementation of SNR-based merging  
- Support further research in multi-source data fusion  

---

## Citation

If you use this code, please cite:

Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022)  
Rethinking Satellite Data Merging: From Averaging to SNR Optimization  
IEEE Transactions on Geoscience and Remote Sensing  
https://doi.org/10.1109/TGRS.2021.3107028

---

## Author

Seokhyeon Kim  
Kyung Hee University  
Water Resources and Environment Big Data Lab (WreBigDL)
