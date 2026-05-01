dataGEN <- function(n, p, ecc, SNRdB) {
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
  
  # Estimation
  SNR <- 10^(SNRdB/10)
  b <- combn(1:p, 2)
  
  # Error covariance matrix: EeeT
  EeeT <- diag(runif(p)) # error variances
  for (i in 1:ncol(b)) {
    EeeT[b[1,i], b[2,i]] <- sqrt(EeeT[b[1,i], b[1,i]] * EeeT[b[2,i], b[2,i]]) * ecc
  }
  
  # Flip the matrix to make it symmetric
  EeeT <- (EeeT + t(EeeT)) - diag(diag(EeeT))
  
  # Generating y and e
  Ey2 <- mean(diag(EeeT) * SNR) # Signal power based on given SNR and EeeT
  
  m <- matrix(0, p + 1, p + 1) # covariance matrix of [y, e]
  m[1,1] <- Ey2
  m[2:(p + 1), 2:(p + 1)] <- EeeT
  
  ye <- matrix(rnorm(n * (p + 1)), n, p + 1) # [y, e]
  ye <- scale(ye, center = TRUE)
  ye <- ye %*% chol(cov(ye))  # scale by the Cholesky factor of the covariance matrix
  ye <- ye %*% chol(m)  # scale by the Cholesky factor of m
  
  y <- ye[, 1]
  e <- ye[, 2:(p + 1)]
  
  return(list(y = y, e = e))
}


WA <- function(EeeT) {
  # WA function
  # This function is to estimate weight for Weighted Average
  #
  # INPUT
  #   EeeT = error covariance matrix (pxp)
  #   
  # OUTPUT 
  #   u = merging weight (px1)
  eta <- matrix(c(rep(1, ncol(EeeT))), ncol = 1)
  u <- as.numeric(solve(t(eta)%*%solve(EeeT)%*%eta))*(solve(EeeT)%*%eta)
  return(u)
}


SNRopt <- function(N, a) {
  # SNRopt function
  # This function is to estimate weight for SNRopt
  #
  # INPUT
  #   N = noise-to-signal ratio matrix (pxp)
  #   a = scaling factor (px1)
  #
  # OUTPUT 
  #   u = merging weight (px1)
  u <- solve(N + a %*% t(a))%*%a
  
  return(u)
}

maxR <- function(a, ExxT) {
  # maxR function
  # This function is to estimate weight for maximizing Pearson R
  #
  # INPUT
  #   a = scaling factor (px1), equivalent to use a = theta; a = rho.*std(x)'
  #
  # OUTPUT 
  #   u: merging weight (px1)
  
  # Generalized Rayleigh quotient
  A <- a %*% t(a)
  B <- ExxT
  
  # Perform the generalized eigenvalue decomposition
  eigen_result <- geigen(A, B)
  
  # Extract the eigenvectors (stored as columns)
  V <- eigen_result$vectors
  # Get the last (leading) eigenvector
  u <- V[, ncol(V)]
  # Normalize u so that its sum is equal to 1
  u <- matrix(u / sum(u), nrow = nrow(a), ncol = 1, byrow = TRUE)
  
  return(u)
}

ECVest <- function(ExxT) {
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
  
  # Parameters
  p <- nrow(ExxT) # number of products
  beta <- 0.5 * min(diag(ExxT)) # tuning parameter for initial a (should be < any of ExxT diagonals)
  lamda <- 0.01 # learning rate
  iters <- 2000 # number of iterations
  
  # Initialization
  eig_result <- eigen(ExxT - beta * diag(p))
  V <- eig_result$vectors # V: eigen (column) vectors
  D <- eig_result$values # D: eigen values (diagonals)
  theta_est <- V[, 1] * sqrt(D[1]) # initialize theta
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
    theta_est <- theta_est - (lamda /    norm(theta_init, type = "2")) * grad
    # project
    theta_est <- theta_est - sign(theta_est) * sqrt(pmax((theta_est^2 - diag(ExxT)), 0))
  }
  
  EeeT_est <- ExxT - theta_est %*% t(theta_est)
  rho2_est <- (theta_est^2) / diag(ExxT)
  
  return(list(EeeT_est = EeeT_est, theta_est = theta_est, rho2_est = rho2_est))
}


SNRest <- function(ExxT, Ey2) {
  # SNRest function
  # This function is for estimating N and a
  #
  # INPUT
  #   ExxT = covariance matrix of x (pxp)
  #   Ey2 = signal power (scalar)
  #
  # OUTPUT
  #   N_est = estimated noise-to-signal ratio (pxp)
  #   a_est = estimated scaling factor (px1)
  
  # Parameters
  C <- ExxT / Ey2
  P <- nrow(C)
  beta <- 0.5 * min(diag(C))
  lambda <- 0.01
  iters <- 2000
  
  # Initialization
  eig_result <- eigen(C - beta * diag(P))
  V <- eig_result$vectors
  D <- eig_result$values
  a_est <- V[, 1] * sqrt(D[1])
  a_est <- a_est * sign(a_est)
  a_init <- a_est
  
  # Iterations
  for (i in 1:iters) {
    grad <- rep(0, P)
    for (j in 1:P) {
      for (k in 1:P) {
        grad[j] <- grad[j] + (j != k) * a_est[k] * sign(a_est[j] * a_est[k] - C[j, k])
      }
    }
    # Descent
    a_est <- a_est - (lambda / norm(a_init, type = "2")) * grad
    
    # Project
    a_est <- a_est - sign(a_est) * sqrt(pmax((a_est ^ 2 - diag(C)), 0))
  }
  
  # Calculate N_est
  N_est <- C - a_est %*% t(a_est)
  a_est <- matrix(a_est, nrow = P, ncol = 1, byrow = TRUE)
  
  return(list(N_est = N_est, a_est = a_est))
}


# EeeTGEN function
# This function is to error covarance matrix
#
# INPUT
#   p = number of datasets (scalar)
#   ecc = error cross-correlation (scalar, [0,1])
#
# OUTPUT
#   EeeT = Error covarance matrix (pxp)
EeeTGEN <- function(p, ecc) {
  b <- combn(p, 2) # combination pairs, e.g., for p=3, (0,1), (0,2), (1,2)
  
  # Error covariance matrix: EeeT
  EeeT <- diag(runif(p)) # error variances
  for (i in 1:ncol(b)) {
    EeeT[b[1,i], b[2,i]] <- sqrt(EeeT[b[1,i], b[1,i]] * EeeT[b[2,i], b[2,i]]) * ecc
  }
  
  EeeT <- (EeeT + t(EeeT)) - diag(diag(EeeT)) # flipping
  
  return(EeeT)
}
