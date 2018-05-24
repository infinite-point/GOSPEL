graph_guided_fused_lasso <- function(X, Y, R, lambda_1, lambda_2,
C, CNorm, maxiter = 1e2, tol = 1e-3, b0, mu = 1e-1){
# Graph guided fused lasso
# Smoothing Proximal Gradient based on FISTA for medium-scale problem with a pre-computed Lipschitz constant
# Reference:
# SMOOTHING PROXIMAL GRADIENT METHOD FOR GENERAL STRUCTURED SPARSE REGRESSION
# Xi Chen, Qihang Lin, Seyoung Kim, Jaime G. Carbonell and Eric P. Xing
# The Annals of Applied Statistics, 2012 Vol. 6, No. 2, 719-752.
# This R code was written by Teppei Shimamura.
# Date: 2017/07/14
# X: inputs
# Y: output
# R: correlation matrix
# lambda_1:  regularization for group norm
# lambda_2: regularization for L1 penalty 
# C:  \sum_|g| by J matrix or |E| by p matrix
# CNorm: ||C||^2
# maxiter: maximum number of iterations
# tol: tau for stopping iterations
# b0: initial coefficients

  n <- nrow(X)
  p <- ncol(X)
  if(!missing(b0)){
    beta <- b0
  } else {
    beta <- Matrix(0,p,1,sparse=TRUE)
  }

  if(missing(C) | missing(CNorm)){
    junk <- generate_C(X, R)
    C <- junk$C
    CNorm <- junk$CNorm
  }

  obj <- rep(0, maxiter)
  density <- rep(0, maxiter)

  C <- C * lambda_1
  w <- beta
  theta <- 1

  XX <- t(X) %*% X
  XY <- t(X) %*% Y

  if(p < 10000) {
    L <- eigs(XX, 1)$values + lambda_1^2 * CNorm / mu
  } else {
    L <- sum(XX^2) + lambda_1^2 * CNorm / mu
  }

  start_time<- proc.time()[3]

  for(iter in 1:maxiter){
    # Update alpha
    A <- hard_threshold(as.matrix(C %*% w / mu), 1)
    # Update gradient
    if(p < 2*n && p < 10000){
      grad <- XX %*% w - XY + t(C) %*% A
    } else {
      grad <- t(X) %*% (X %*% w - Y) + t(C) %*% A
    }
    beta_new <- soft_threshold(w - (1/L) * grad, lambda_2 / L)
    # Update beta
    density[iter] <- sum(beta_new != 0) / p
    # Update theta
    theta_new <- (sqrt(theta^4 + 4 * theta^2) - theta^2) / 2
    # Update w
    w <- beta_new + (1 - theta) / theta * theta_new * (beta_new - beta)
    # Calculate object function
    obj[iter] <- sum((Y - X %*% beta_new)^2) / 2 + sum(abs(C %*% beta_new)) + lambda_2 * sum(abs(beta_new))
    # Check convergence
    if(iter > 10 && (abs(obj[iter] - obj[iter - 1]) / abs(obj[iter - 1]) < tol)) break
    beta <- beta_new
    theta <- theta_new
  }

  end_time<- proc.time()[3]

  return(list(beta = beta, theta = theta, obj = obj[1:iter],
                 density = density[1:iter], iter = iter,
                 time = end_time - start_time))

}