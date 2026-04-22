# SNRopt: Satellite Data Merging using SNR Optimization

This repository provides a modular Python implementation of data merging methods based on:

Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021).  
Rethinking Satellite Data Merging: From Averaging to SNR Optimization.  
IEEE Transactions on Geoscience and Remote Sensing

The code includes implementations of:
- Weighted Average (WA)
- SNR-optimal merging (SNRopt)
- Maximum Correlation merging (maxR)
- Noise-Constrained (NC) estimation
- SNR-based parameter estimation (SNRest)

---

## Files

- `SNRopt_modular.ipynb`  
  Example workflow (main notebook)

- `snropt_helpers.py`  
  Core functions (modular implementation)

---

## Requirements

```bash
pip install numpy matplotlib
