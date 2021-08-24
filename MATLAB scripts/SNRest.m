function [N_est,a_est] = SNRest(ExxT, Ey2)
%% This function is for eatimating N and a
%
% INPUT
%   ExxT = covariance matrix of x (pxp)
%   Ey2 = signal power (scalar)
%
% OUTPUT
%   N_est = estimated noise-to-signal ratio (pxp)
%   a_est = estimated scaling factor (px1)
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
% Parameters
C = ExxT/Ey2;
P = size(C,1); % number of products
beta = 0.5*min(diag(C)); % tuning parameter for itial a (should be < any of C diagonals)
% beta = 0.1; % tuning parameter for unitializing a (should be < any of C diagonals)
lamda = 0.01; % learning rate
iters = 2000; % number of iterations

% Initialization
[V,D] = eig(C-beta*eye(P));% V: eigen (column) vectors; D: eigen values (diagonals)
a_est = V(:,end)*sqrt(D(end)); % initialize a
a_est = a_est .* sign(a_est);  % assuming '+'sign'
a_init = a_est;

% Iterations
for i = 1:iters
    grad = zeros(P,1);
    for j = 1:P
        for k = 1:P
            grad(j) = grad(j) + (j~=k) * a_est(k)*sign(a_est(j)*a_est(k) - C(j,k));
        end
    end
    % descent
    a_est = a_est - (lamda/norm(a_init))*grad;
    % project
    a_est = a_est - sign(a_est) .* sqrt(max((a_est.^2 - diag(C)),0));
    
end

N_est = C-a_est*a_est';

end