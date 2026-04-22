function u = maxR(a, ExxT)
% This function is to estimate weight for maximizing Pearson R
%
% INPUT
%   a = scaling factor (px1), equivalent to use a = theta; a = rho.*std(x)'
%
% OUTPUT 
%   u: merging weight (px1)
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
% Generalized Rayleigh quotient
% A*u = lamda*B*u
A = a*a';
B = ExxT;

[V,~] = eig(A,B);% V: eigen (column) vectors
u = V(:,end); % leading eigen vector
u = u./sum(u); % normalizing sum to 1

end