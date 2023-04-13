function [EeeT_est, theta_est, rho2_est] = ECVest(ExxT)
%% This function is a modified version of SNRest to estimate TC-like results
% Version 2 (updated on 13 Apr 2023)
% Update details: same results but with more detailed steps
%
% INPUT
%   ExxT = covariance matrix of x (p x p)
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
p = size(ExxT,1); % number of products
beta = 0.5*min(diag(ExxT)); % tuning parameter for initial a (should be < any of ExxT diagonals)
% beta = 0.1; % tuning parameter for uninitializing a (should be < any of ExxT diagonals)
lamda = 0.01; % learning rate
iters = 2000; % number of iterations

% Initialization
[V,D] = eig(ExxT-beta*eye(p)); % V: eigen (column) vectors; D: eigen values (diagonals)
theta_est = V(:,end)*sqrt(D(end)); % initialize theta
theta_est = theta_est .* sign(theta_est);  % assuming '+' sign
theta_init = theta_est;

% Iteration
for i = 1:iters
    % Calculate gradient of the cost function
    grad = zeros(p,1);
    for j = 1:p
        % Calculate the inner summation term of the gradient for j-th element
        grad_j = 0;
        for k = 1:p % Loop through the indices of the gradient vector again            
            if j~=k % Only update the gradient if j is not equal to k (off-diagonals)
                prod = theta_est(j)*theta_est(k); % Calculate the product of the j-th and k-th entries of theta_est
                diff = prod - ExxT(j,k); % Calculate the difference between the product and the corresponding entry of ExxT                
                grad_j = grad_j + theta_est(k)*sign(diff); % Add the contribution of the k-th entry to the j-th element of the gradient
            end
        end        
        grad(j) = grad_j; % Store the value of the j-th element of the gradient
    end
    % Update theta_est using gradient descent
    theta_est = theta_est - (lamda/norm(theta_init))*grad;

    % Projecting non-negativity constraint on the updated theta_est values
    theta_est = theta_est - sign(theta_est) .* sqrt(max((theta_est.^2 - diag(ExxT)),0));         
    
end

% Calculate EeeT_est and rho2_est using theta_est
EeeT_est = ExxT - theta_est*theta_est';
rho2_est = (theta_est.^2)./diag(ExxT);

end