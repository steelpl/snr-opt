function [y,e] = dataGEN(n,p,ecc,SNRdB)
% This function is to generate orthogonal y and e
%
% INPUT
%   n = data length (scalar)
%   p: number of datasets (scalar)
%   ecc = error cross-correlation (scalar, [0,1])
%   SNRdB = signal-to-noise ratio in dB (scalar)
%
% OUTPUT
%   y = signal (nx1) zero mean
%   e = error (nxp) zero mean 
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
SNR = 10^(SNRdB/10);
b = nchoosek(1:p,2); % combination pairs e.g. for p=3, (1,2), (1,3), (2,3)

% Error covarance matrix: EeeT
EeeT = diag(rand(p,1)); % error variances
for i = 1:size(b,1)
    EeeT(b(i,1),b(i,2)) = sqrt(EeeT(b(i,1),b(i,1))*EeeT(b(i,2),b(i,2)))*ecc;
    
end
EeeT = (EeeT+EeeT') - eye(p).*diag(EeeT); % flipping

% Generating y and e
Ey2 = mean(diag(EeeT)*SNR); % Signal power based on given SNR and EeeT

m = zeros(p+1); % covariance matrix of [y, e]
m(1,1) = Ey2; m(2:end,2:end) = EeeT;

ye = randn(n, p+1); %[y, e]
ye = ye - mean(ye); % mean removal
ye = ye/chol(cov(ye));
ye = ye * chol(m);

y = ye(:,1);
e = ye(:,2:end);

end