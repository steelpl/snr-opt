function u = WA(EeeT)
%% This function is to estimate weight for Weighted Average
%
% INPUT
%   EeeT = error covariance matrix (pxp)
%   
% OUTPUT 
%   u = merging weight (px1)
%
% REFERENCE
% For more details, see:
%
% Kim, S., Sharma, A., Liu, Y. Y., & Young, S. I. (2021). 
% Rethinking Satellite Data Merging: From Averaging to SNR Optimization.
% IEEE Trans Geosci Remote Sens
%
% If you use the methods presented in the paper and/or this example, 
% please cite this paper where appropriate.
%
%% Estimation
eta = ones(size(EeeT,2),1);
u = (eta'*(EeeT\eta))\(EeeT\eta);

end