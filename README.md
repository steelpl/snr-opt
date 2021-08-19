# SNR-opt

This is for sharing codes (MATLAB) used in 
> Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021). "Rethinking Satellite Data Merging: From Averaging to SNR Optimization." IEEE
Trans Geosci Remote Sens, Accepted

The final version will be available soon on the journal website, and a preprint that is not much different from the final version can be obtained from [ResearchGate](https://www.researchgate.net/publication/349961492_Rethinking_Satellite_Data_Merging_From_Averaging_to_SNR_Optimization_techrxiv14214035) or [Techrxiv](https://www.techrxiv.org/articles/preprint/Rethinking_Satellite_Data_Merging_From_Averaging_to_SNR_Optimization/14214035)

## Code description
The codes consist of **nine MATLAB** scripts and a brief explanation for each one is as follows.
### Two main scripts in which four merging methods (weighted averaging, SNR-opt, maxR, unweighted averaging) are compared
1. *mergingExample_the.m*: an example to merge three (or any number) synthetic zero mean datasets (no time series).
2. *mergingExample_syn.m*: an example to merge three (or any number) synthetic zero mean datasets (time series).

### Seven functions used in the main scripts
1. *WA.m*: to calculate merging coefficients for weighted averaging.
2. *SNRopt.m*: to calculate merging coefficients for SNR-opt.
3. *maxR.m*: to calculate merging coefficients for maxR (maximizing correlation).
4. *SNRest.m*: to estimate noise-to-signal ratio and scaling factor for given covariance matrix and signal power (any number of products).
5. *ECVest.m*: to estimate error covariance matrix and data-truth correlation for given covariance matrix (any number of products). This is equivalent to **Triple Collocation** but produces no failures and is applicable for any number of products.
6. *EeeTGEN.m*: to generate error covariance matrix for given number of products and error cross-correlation
7. *dataGEN.m*: to generate orthogonal y (truth) and e (error) for given data length, number of products, error cross-correlation, and signal-to-noise ratio
  
## Disclaimer
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
