function [N_est,a_est] = SNRest(ExxT, Ey2)
%% This function is for estimating N and a
% Version 2 (updated on 13 Apr 2023)
% Update details: same results but with more detailed steps
%
% INPUT
%   ExxT = covariance matrix of x (p x p)
%   Ey2 = signal power (scalar)
%
% OUTPUT
%   N_est = estimated noise-to-signal ratio (p x p)
%   a_est = estimated scaling factor (p x 1)
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
C = ExxT/Ey2; % Calculate normalized covariance matrix C by dividing ExxT by Ey2
p = size(C,1); % Determine the number of products (size of C matrix)
beta = 0.5*min(diag(C)); % Calculate initial value for beta parameter (tuning parameter for initializing a), which is set to 0.5*min(diag(C))
lamda = 0.01; % Set learning rate to 0.01
iters = 2000; % Set number of iterations for gradient descent to 2000

% Initialization
[V,D] = eig(C-beta*eye(p)); % Compute eigenvectors V and eigenvalues D of the matrix C-beta*I, where I is identity matrix of size P.
a_est = V(:,end)*sqrt(D(end)); % Initialize a_est to the eigenvector corresponding to the largest eigenvalue of C-beta*I, scaled by the square root of that eigenvalue.
a_est = a_est .* sign(a_est); % Make all entries of a_est positive by taking the element-wise sign of a_est
a_init = a_est; % Store initial value of a_est

% Iterations
for i = 1:iters % Loop through the number of iterations specified
    grad = zeros(p,1); % Initialize gradient to zero
    for j = 1:p % Loop through the indices of the gradient vector
        grad_j = 0; % Initialize the j-th element of the gradient to zero
        for k = 1:p % Loop through the indices of the gradient vector again
            if j~=k % Only update the gradient if j is not equal to k (off-diagonals)
                prod = a_est(j)*a_est(k); % Calculate the product of the j-th and k-th entries of a_est
                diff = prod - C(j,k); % Calculate the difference between the product and the corresponding entry of C
                grad_j = grad_j + a_est(k)*sign(diff); % Add the contribution of the k-th entry to the j-th element of the gradient
            end
        end
        grad(j) = grad_j; % Store the value of the j-th element of the gradient
    end
    % Update a_est using gradient descent
    a_est = a_est - (lamda/norm(a_init))*grad; 

    % Projecting non-negativity constraint on the updated a_est values
    a_est = a_est - sign(a_est) .* sqrt(max((a_est.^2 - diag(C)),0));

end

N_est = C-a_est*a_est'; % Calculate the estimated noise-to-signal ratio using the formula in the paper

end