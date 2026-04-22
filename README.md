# # SNR-opt: Reference-Free Multi-Source Data Fusion Beyond TC-Based Merging

This repository provides MATLAB and Python implementations accompanying:

Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).  
**Rethinking satellite data merging: from averaging to SNR optimization.**  
*IEEE Transactions on Geoscience and Remote Sensing*, 60, 1–15.  
https://doi.org/10.1109/TGRS.2021.3107028

---

## Overview

This repository implements and compares multiple data merging methods for multi-source observations:

- Weighted Averaging (WA)
- SNR-optimal merging (SNRopt)
- Correlation-maximizing merging (maxR)
- Equal Weighting (EW)

The implementations include:
- Synthetic data generation with controlled error structure
- Covariance-based parameter estimation
- Performance evaluation using MSE and R²

---

## Signal + Noise Model

$$
x = a y + e
$$

$$
Q = E[xx^T] = a a^T E[y^2] + E[ee^T]
$$

---

## SNR-optimal Merging

$$
u = (N + a a^T)^{-1} a
$$

where:

$$
C = \frac{Q}{E[y^2]} = a a^T + N
$$

---

## Repository Structure

```text
├── Jupyter Notebook/
│ ├── SNRopt_modular.ipynb
│ └── snropt_helpers.py
├── MATLAB scripts/
│ ├── mergingExample_the.m
│ ├── mergingExample_syn.m
│ ├── WA.m
│ ├── SNRopt.m
│ ├── maxR.m
│ ├── SNRest.m
│ ├── NC.m
│ ├── EeeTGEN.m
│ └── dataGEN.m
├── README.md
└── LICENSE
```

## Usage

### MATLAB

See detailed documentation:

➡️ `MATLAB scripts/README.md`

Includes:
- theoretical example (covariance-based)
- synthetic time-series example
- full comparison of WA, SNRopt, maxR, EW

---

### Python

See detailed documentation:

➡️ `Jupyter Notebook/README.md`

Includes:
- modular implementation
- interactive notebook workflow
- reproducible experiments

---

## Key Insights

- SNRopt provides **minimum MSE** when parameters are known
- Performance is **sensitive to estimation errors**
- WA is **more robust but suboptimal**
- NC generalizes TC and avoids failure cases

---

## Purpose

This repository is intended to:

- Reproduce results from the paper
- Provide transparent implementations (MATLAB & Python)
- Support research in multi-source data fusion

---

## Citation

If you use this code, please cite:

Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2022).  
**Rethinking satellite data merging: from averaging to SNR optimization.**  
*IEEE Transactions on Geoscience and Remote Sensing*, 60, 1–15.  
https://doi.org/10.1109/TGRS.2021.3107028

---

## Author

Seokhyeon Kim  
(on behalf of all coauthors)  
Kyung Hee University  
Water Resources and Environment Big Data Lab (WreBigDL)  
https://sites.google.com/view/wrebigdl/
