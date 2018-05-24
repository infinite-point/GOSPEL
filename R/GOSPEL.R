require("data.table")
require("rARPACK")
require("igraph")
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
