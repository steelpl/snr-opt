function [errVar_EIVD, rho2_EIVD,ECC2_EIVD] = EIVD(DS)

% Seokhyeon Kim
% Date: 15 Dec, 2019

% Catch errors in inputs
if size(DS,2) ~= 3
    error('Error: Input data X must be an N x 3 matrix');
end

if any(isnan(DS(:)))
    error('Error: Input data X must not contain NaNs');
end

if length(unique(DS(:,1))) == 1 || length(unique(DS(:,2))) == 1 || length(unique(DS(:,3))) == 1
    error('Error: the sample variance of each of the columns of X must be non-zero. Increase your sample size or reconsider using IVd.');
end

% three input variables - note the sequence
x = DS(:,3); % ds3 -> x
z = DS(:,1); % ds1 -> z
w = DS(:,2); % ds2 -> w

% % Generate double instrumental variables (lag-1 [day] time series)
xLag1 = x; xLag1(1:end-1) = x(2:end);
zLag1 = z; zLag1(1:end-1) = z(2:end);
wLag1 = w; wLag1(1:end-1) = w(2:end);

% Cross- and auto-covariances
Cxx = var(x); Cxx = Cxx(1,1);
Czz = cov(z,z); Czz = Czz(1,1);
Cww = cov(w,w); Cww = Cww(1,1);

Czw = cov(z,w); Czw = Czw(1,2);
Cxz = cov(x,z); Cxz = Cxz(1,2);
Cxw = cov(x,w); Cxw = Cxw(1,2);

Lxx = cov(x,xLag1); Lxx = Lxx(1,2);
Lzz = cov(z,zLag1); Lzz = Lzz(1,2);
Lww = cov(w,wLag1); Lww = Lww(1,2);
Lxx = abs(Lxx); Lzz = abs(Lzz); Lww = abs(Lww); % absolute values

% A matrix
A = [[1,0,0,0,1,0,0,0];
    [0,1,0,0,0,1,0,0]
    [0,0,1,0,0,0,1,0]
    [0,0,0,1,0,0,0,1]
    [1,0,0,0,0,0,0,0]
    [1,0,0,0,0,0,0,0]
    [0,1,0,0,0,0,0,0]
    [0,0,1,0,0,0,0,0]
    [0,0,0,1,0,0,0,0]
    [0,0,0,1,0,0,0,0]];

% b matrix
b = [Cxx;Czz;Cww;Czw;...
    Cxz*sqrt(Lxx/Lzz);...
    Cxw*sqrt(Lxx/Lww);...
    Cxz*sqrt(Lzz/Lxx);...
    Cxw*sqrt(Lww/Lxx);...
    Cxw*sqrt(Lzz/Lxx);...
    Cxz*sqrt(Lww/Lxx)];

% error metrix (e) to be calcualted by OLS
e = (A'*A)\A'*b;

errVar_x = e(5,1);
errVar_z = e(6,1);
errVar_w = e(7,1);
errVar_EIVD = [errVar_x;errVar_z;errVar_w];

rho2_x = e(1,1)/(e(1,1)+e(5,1));
rho2_z = e(2,1)/(e(2,1)+e(6,1));
rho2_w = e(3,1)/(e(3,1)+e(7,1));
rho2_EIVD = [rho2_x;rho2_z;rho2_w];

ECC2_zw = e(8,1)^2/(e(6,1)*e(7,1));
ECC2_EIVD = [0;0;ECC2_zw];

% If errVars are negative, display a warning.
if any(errVar_EIVD < 0)
    warning('Warning: at least one calculated errVar is negative. This can happen if the sample size is too small, or if one of the assumptions of IVd is violated.');
end
% If squared correlation coefficients are negative, display a warning.
if any(rho2_EIVD < 0)
    warning('Warning: at least one calculated squared correlation coefficient is negative. This can happen if the sample size is too small, or if one of the assumptions of IVd is violated.');
end
% If squared correlation coefficients are negative, display a warning.
if any(ECC2_zw < 0)
    warning('Warning: squared ECC_zw is negative. This can happen if the sample size is too small, or if one of the assumptions of IVd is violated.');
end

end