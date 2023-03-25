# dataGEN function
# This function is to generate orthogonal y and e
#
# INPUT
#   n = data length (scalar)
#   p: number of datasets (scalar)
#   ecc = error cross-correlation (scalar, [0,1])
#   SNRdB = signal-to-noise ratio in dB (scalar)
#
# OUTPUT
#   y = signal (nx1)
#   e = error (nxp)
dataGEN <- function(n, p, ecc, SNRdB) {
  # Estimation
  SNR <- 10^(SNRdB/10)
  b <- t(which(upper.tri(matrix(1, p, p), diag = FALSE))) # combination pairs
  
  # Error covariance matrix: EeeT
  EeeT <- diag(runif(p)) # error variances
  for (i in 1:nrow(b)) {
    EeeT[b[i,1], b[i,2]] <- sqrt(EeeT[b[i,1], b[i,1]] * EeeT[b[i,2], b[i,2]]) * ecc
    EeeT[b[i,2], b[i,1]] <- EeeT[b[i,1], b[i,2]]
  }
  
  # Generating y and e
  Ey2 <- mean(diag(EeeT) * SNR) # Signal power based on given SNR and EeeT
  
  m <- matrix(0, p + 1, p + 1) # covariance matrix of [y, e]
  m[1,1] <- Ey2
  m[2:(p + 1), 2:(p + 1)] <- EeeT
  
  ye <- matrix(rnorm(n * (p + 1)), n, p + 1) # [y, e]
  ye <- ye - apply(ye, 2, mean) # mean removal
  cov <- cov(t(ye)) # compute covariance matrix of transposed array
  cholesky_cov <- chol(cov) # Cholesky decomposition of covariance matrix
  ye <- ye %*% solve(cholesky_cov) # normalize by inverse of Cholesky decomposition of covariance matrix
  ye <- ye %*% chol(m) # multiply by Cholesky decomposition of m
  
  y <- ye[, 1]
  e <- ye[, 2:(p + 1)]
  
  return(list(y = y, e = e))
}

# WA function
# This function is to estimate weight for Weighted Average
#
# INPUT
#   EeeT = error covariance matrix (pxp)
#   
# OUTPUT 
#   u = merging weight (px1)
WA <- function(EeeT) {
  eta <- matrix(1, nrow = nrow(EeeT), ncol = 1)
  u <- solve(EeeT) %*% eta / t(eta) %*% solve(EeeT) %*% eta
  
  return(u)
}

# SNRopt function
# This function is to estimate weight for SNRopt
#
# INPUT
#   N = noise-to-signal ratio matrix (pxp)
#   a = scaling factor (px1)
#
# OUTPUT 
#   u = merging weight (px1)
SNRopt <- function(N, a) {
  u <- solve(N + outer(a, a)) %*% a
  
  return(u)
}

# maxR function
# This function is to estimate weight for maximizing Pearson R
#
# INPUT
#   a = scaling factor (px1), equivalent to use a = theta; a = rho.*std(x)'
#
# OUTPUT 
#   u: merging weight (px1)
maxR <- function(a, ExxT) {
  # Generalized Rayleigh quotient
  # A * u = lamda * B * u
  A <- outer(a, a)
  B <- ExxT
  
  eig <- eigen(solve(B) %*% A)
  u <- eig$vectors[, which.max(eig$values)]
  u <- u / sum(u) # normalizing sum to 1
  
  return(u)
}

# ECVest function
# This function is a modified version of SNRest to estimate TC-like results
#
# INPUT
#   ExxT = covariance matrix of x (pxp)
#
# OUTPUT 
#   EeeT_est = estimated error covariance matrix (pxp)
#   theta_est = estimated theta (px1, = a*sqrt(Ey2))
#   rho2_est = estimated squared data-truth correlation (px1)
ECVest <- function(ExxT) {
  # Parameters
  p <- nrow(ExxT) # number of products
  beta <- 0.5 * min(diag(ExxT)) # tuning parameter for itial a (should be < any of ExxT diagonals)
  # beta <- 0.1 # tuning parameter for unitializing a (should be < any of ExxT diagonals)
  lamda <- 0.01 # learning rate
  iters <- 2000 # number of iterations
  
  # Initialization
  eig <- eigen
  V <- eig$vectors
  D <- diag(eig$values)
  max_idx <- which.max(eig$values)
  theta_est <- V[, max_idx] * sqrt(D[max_idx]) # initialize theta
  theta_est <- theta_est * sign(theta_est) # assuming '+' sign'
  theta_init <- theta_est
  
  # Iteration
  for (i in 1:iters) {
    grad <- rep(0, p)
    for (j in 1:p) {
      for (k in 1:p) {
        grad[j] <- grad[j] + (j != k) * theta_est[k] * sign(theta_est[j] * theta_est[k] - ExxT[j, k])
      }
    }
    # descent
    theta_est <- theta_est - (lamda / norm(theta_init)) * grad
    # project
    theta_est <- theta_est - sign(theta_est) * sqrt(pmax(theta_est^2 - diag(ExxT), 0))
  }
  
  EeeT_est <- ExxT - outer(theta_est, theta_est)
  rho2_est <- theta_est^2 / diag(ExxT)
  
  return(list(EeeT_est = EeeT_est, theta_est = theta_est, rho2_est = rho2_est))
}

# SNRest function
# This function is for eatimating N and a
#
# INPUT
#   ExxT = covariance matrix of x (pxp)
#   Ey2 = signal power (scalar)
#
# OUTPUT
#   N_est = estimated noise-to-signal ratio (pxp)
#   a_est = estimated scaling factor (px1)
SNRest <- function(ExxT, Ey2) {
  # Parameters
  C <- ExxT / Ey2
  P <- nrow(C) # number of products
  beta <- 0.5 * min(diag(C)) # tuning parameter for itial a (should be < any of C diagonals)
  # beta <- 0.1 # tuning parameter for unitializing a (should be < any of C diagonals)
  lamda <- 0.01 # learning rate
  iters <- 2000 # number of iterations
  
  # Initialization
  eig <- eigen(C - beta * diag(P))
  a_est <- eig$vectors[, P] * sqrt(eig$values[P]) # initialize a
  a_est <- a_est * sign(a_est) # assuming '+' sign'
  a_init <- a_est
  
  # Iterations
  for (i in 1:iters) {
    grad <- rep(0, P)
    for (j in 1:P) {
      for (k in 1:P) {
        grad[j] <- grad[j] + (j != k) * a_est[k] * sign(a_est[j] * a_est[k] - C[j, k])
      }
    }
    # descent
    a_est <- a_est - (lamda / norm(a_init)) * grad
    # project
    a_est <- a_est - sign(a_est) * sqrt(pmax(a_est^2 - diag(C), 0))
  }
  
  N_est <- C - outer(a_est, a_est)
  
  return(list(N_est = N_est, a_est = a_est))
}