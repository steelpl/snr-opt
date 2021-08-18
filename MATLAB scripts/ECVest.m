function [EeeT_est,theta_est,rho2_est] = ECVest(ExxT)
%% This function is a modified version of SNRest to estimate TC-like results
%
% INPUT
%   ExxT = covariance matrix of x (pxp)
%
% OUTPUT 
%   EeeT_est = estimated error covariance matrix (pxp)
%   theta_est = estimated theta (px1, = a*sqrt(Ey2))
%   rho2_est = estimated squared data-truth correlation (px1)
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
p = size(ExxT,1); % number of products
beta = 0.5*min(diag(ExxT)); % tuning parameter for itial a (should be < any of ExxT diagonals)
% beta = 0.1; % tuning parameter for unitializing a (should be < any of ExxT diagonals)
lamda = 0.01; % learning rate
iters = 2000; % number of iterations

% Initialization
[V,D] = eig(ExxT-beta*eye(p)); % V: eigen (column) vectors; D: eigen values (diagonals)
theta_est = V(:,end)*sqrt(D(end)); % initialize theta
theta_est = theta_est .* sign(theta_est);  % assuming '+'sign'
theta_init = theta_est;

% Iteration
for i = 1:iters
    grad = zeros(p,1);
    for j = 1:p
        for k = 1:p
            grad(j) = grad(j) + (j~=k) * theta_est(k)*sign(theta_est(j)*theta_est(k) - ExxT(j,k));
        end
    end
    % descent
    theta_est = theta_est - (lamda/norm(theta_init))*grad;
    % project
    theta_est = theta_est - sign(theta_est) .* sqrt(max((theta_est.^2 - diag(ExxT)),0));
    
end
EeeT_est = ExxT-theta_est*theta_est';
rho2_est = (theta_est.^2)./diag(ExxT);

end