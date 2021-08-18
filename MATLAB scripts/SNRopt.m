function u = SNRopt(N,a)
%% This function is to estimate weight for SNRopt
%
% INPUT
%   N = noise-to-signal ratio matrix (pxp)
%   a = scaling factor (px1)
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
u = (N+a*a')\a;

end

