function EeeT = EeeTGEN(p,ecc)
% This function is to error covarance matrix
%
% INPUT
%   p = number of datasets (scalar)
%   ecc = error cross-correlation (scalar, [0,1])
%
% OUTPUT
%   EeeT = Error covarance matrix (pxp)
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
b = nchoosek(1:p,2); % combination pairs e.g. for p=3, (1,2), (1,3), (2,3)

% Error covarance matrix: EeeT
EeeT = diag(rand(p,1)); % error variances
for i = 1:size(b,1)
    EeeT(b(i,1),b(i,2)) = sqrt(EeeT(b(i,1),b(i,1))*EeeT(b(i,2),b(i,2)))*ecc;
    
end
EeeT = (EeeT+EeeT') - eye(p).*diag(EeeT); % flipping

end