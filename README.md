# SNR-opt

This is for sharing codes (MATLAB) used in 
> S. Kim, A. Sharma, Y. Y. Liu and S. I. Young, "Rethinking Satellite Data Merging: From Averaging to SNR Optimization," in IEEE Transactions on Geoscience and Remote Sensing, vol. 60, pp. 1-15, 2022, Art no. 4405215, doi: http://dx.doi.org/10.1109/TGRS.2021.3107028

## Code description
The codes consist of **nine MATLAB** scripts and a brief explanation for each one is as follows.
### Two main scripts in which four merging methods (weighted averaging, SNR-opt, maxR, unweighted averaging) are compared
1. ***mergingExample_the.m***: an example to merge three (or any number) synthetic zero mean datasets (no time series).
2. ***mergingExample_syn.m***: an example to merge three (or any number) synthetic zero mean datasets (time series).

### Seven functions used in the main scripts
1. ***WA.m***: to calculate merging coefficients for weighted averaging.
2. ***SNRopt.m***: to calculate merging coefficients for SNR-opt.
3. ***maxR.m***: to calculate merging coefficients for maxR (maximizing correlation).
4. ***SNRest.m***: to estimate noise-to-signal ratio and scaling factor for given covariance matrix and signal power (any number of products).
5. ***ECVest.m***: to estimate error covariance matrix and data-truth correlation for given covariance matrix (any number of products). This is equivalent to **Triple Collocation** but produces no failures and is applicable for any number of products.
6. ***EeeTGEN.m***: to generate error covariance matrix for given number of products and error cross-correlation
7. ***dataGEN.m***: to generate orthogonal y (truth) and e (error) for given data length, number of products, error cross-correlation, and signal-to-noise ratio
