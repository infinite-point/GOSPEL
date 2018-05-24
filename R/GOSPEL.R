require("rARPACK")
require("bayesopt")
require("Matrix")

GOSPEL <- function(Ker_X, Ker_Y, CKA_mat, lambda_1 = 0.01, tau, lambda_2){
  if(missing(tau)) tau <- quantile(unique(CKA_mat[CKA_mat!=1]),probs=seq(from=0.75,to=0.95,length=10))
  if(missing(lambda_2)) lambda_2 <- 10^(seq(from=-2,to=2,length=10))
  objective_func <- function(x,y,z){
    n <- nrow(CKA_mat)
    R <- matrix(0,n,n)
    R[CKA_mat >= x] <- 1
   try({
      res <- try(graph_guided_fused_lasso(Ker_X, Ker_Y, R, y, z),silent=FALSE)
      if(inherits(res,"try-error")){ # If graph guided fused lasso cannot be fit, just return Inf
        bic <- Inf
      } else {
        bic <- BIC_GP(Ker_X, Ker_Y, res$beta)
      }
   })
   return(-bic)
  }
  res1 <- try(bayesopt(objective_func,tau,lambda_2,lambda_1,kernel="square_exp",init_size=10),silent=FALSE)
  if(all(res1$ys==-Inf)){ # If bayesopt fails, return NA
    result <- list(beta = rep(NA,ncol(Ker_X)),
                   bic = NA,
                   probs = NA,
                   lambda_2 = NA,
                   lambda_1 = NA)
    return(result)
  }
  n <- nrow(CKA_mat)
  R <- matrix(0,n,n)
  R[CKA_mat >= res1$opt_x[,1]] <- 1
  res2 <- graph_guided_fused_lasso(Ker_X, Ker_Y, R, res1$opt_x[,2], res1$opt_x[,3])
  result <- list(beta = res2$beta, bic = - res1$opt_y, tau = res1$opt_x[,1],
                      lambda_2=res1$opt_x[,2], lambda_1=res1$opt_x[,3], R = R)
  return(result)
}

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

soft_threshold <- function(vv, lambda_2){
  n <- nrow(vv)
  p.ind <- which(vv > lambda_2)
  n.ind <- which(vv < - lambda_2)
  res <- Matrix(0, n, 1, sparse=TRUE)
  res[p.ind] <- vv[p.ind] - lambda_2
  res[n.ind] <- vv[n.ind] + lambda_2
  return(res)
}

hard_threshold <- function(vv, lambda_2){
  n <- nrow(vv)
  p.ind <- which(vv > lambda_2)
  n.ind <- which(vv < - lambda_2)
  res <- Matrix(0, n, 1, sparse=TRUE)
  res[p.ind] <- lambda_2
  res[n.ind] <- - lambda_2
  return(res)
}

gen_kernel_matrix <- function(X, Y, gamma){
  size <- dim(X)
  n <- size[1]
  p <- size[3]
  H <- diag(n) - (1/n) * matrix(1, n, n)
  Ker_X <- array(0, dim = c(n, n, p)) 
  for (i in 1:p){
    XX <- X[, , i] 
    DX <- diag(rowSums(XX))
    LX <- DX - XX   
    Ker_XX <- exp(-LX/gamma)
    Ker_XX <- H %*% Ker_XX %*% H
    norm <- sqrt(sum(Ker_XX^2))
    Ker_X[, , i] <- Ker_XX /norm
    DY <- diag(rowSums(Y))
    LY <- DY - Y     
    Ker_YY <- exp(-LY/gamma)
    Ker_YY <- H %*% Ker_YY %*% H
    norm <- sqrt(sum(Ker_YY^2))
    Ker_Y <- Ker_YY /norm
  }
    return(list(Ker_X = Ker_X, Ker_Y = Ker_Y))
}

calculate_R <- function(X, Y, gamma, verbose = TRUE){
  size <- dim(X)
  n_node <- size[1]
  p <- size[3]
  KK <- gen_kernel_matrix(X, Y, gamma)
  Ker_X <- KK$Ker_X   
  Ker_Y <- KK$Ker_Y 
  CKA_mat <- diag(p)
  for (i in 1:(p-1)){  
    Ker_X1 <- Ker_X[, , i]
    for (j in (i+1):p){
      if(verbose) cat(paste0("Generating a correlation matrix (" ,i,",",j, ")\n"))
      Ker_X2 <- Ker_X[, , j]
      CKA_mat[i, j] <- CKA_mat[j, i] <- abs(sum(Ker_X1 * Ker_X2))
    }
  }
  dim(Ker_X) <- c((dim(Ker_X)[1])^2, p)
  dim(Ker_Y) <- c((dim(Ker_Y)[1])^2, 1)
  return(list(CKA_mat = CKA_mat, Ker_X = Ker_X, Ker_Y = Ker_Y))
}

generate_C <- function(X, R){
# C: edge by vertex incident matrix as in the paper of Chen et al.
# CNorm in Proposition 2 in the paper of Chen et al.
  R[lower.tri(R)] <- 0
  E <- which(R != 0, arr.ind = TRUE)
  nE <- nrow(E)
  p <- ncol(X)
  C <- Matrix(0, nE, p, sparse=TRUE)
  for(i in 1:nE){
    C[i, E[i, 1]] <- C[i, E[i, 2]] <- 1
  }
  CNorm <- 2 * max(rowSums(C^2))
  return(list(C = C, CNorm = CNorm))
}


BIC_GP <- function(Ker_X, Ker_Y, beta){
  n2 <- dim(Ker_X)[1]
  res2 <- mean((Ker_Y - Ker_X %*% beta)^2)
  loglik <- -(n2/2) * (log(2 * pi * res2) + 1)
  df <- length(table(as.vector(beta))) - 1
  bic <- -2 * loglik + log(n2) * df
  return(bic)  
}

