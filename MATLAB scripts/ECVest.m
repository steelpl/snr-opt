function [EeeT_est, theta_est, rho2_est] = ECVest(covx)
%% This function is a modified version of SNRest to estimate TC-like results
% Version 3 (updated on 23 Mar 2024)
% Update details: same results but with more detailed steps
%
% INPUT
%   covx = covariance matrix of x (p x p)
%
% OUTPUT 
%   EeeT_est = estimated error covariance matrix (p x p)
%   theta_est = estimated theta (p x 1, = a*sqrt(Ey2))
%   rho2_est = estimated squared data-truth correlation (p x 1)
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
p = size(covx,1); % number of products
beta = 0.5*min(diag(covx)); % tuning parameter for initial theta (should be < any of covx diagonals)
lamda = 0.01; % learning rate
iters = 2000; % number of iterations

% Initialization
[V,D] = eig(covx-beta*eye(p)); % V: eigen (column) vectors; D: eigen values (diagonals)
theta_est = V(:,end)*sqrt(D(end)); % initialize theta
theta_est = theta_est .* sign(theta_est);  % assuming '+' sign
theta_init = theta_est;
lamda = lamda/norm(theta_init); % normalizing lamda considering the norm of theta_init
eta = ones(p,1);

% Iteration
for i = 1:iters
    % Calculate gradient of the cost function
    grad = ((eta*eta'-eye(p)).*sign(theta_est*theta_est'-covx))*theta_est;

    % Update theta_est using gradient descent
    theta_est = theta_est - lamda*grad;

    % Projecting non-negativity constraint on the updated theta_est values
    theta_est = theta_est - sign(theta_est) .* sqrt(max(diag(theta_est*theta_est'-covx),0));    
    
end

% Calculate EeeT_est and rho2_est using theta_est
EeeT_est = covx - theta_est*theta_est';
rho2_est = diag((theta_est*theta_est')./covx);

end