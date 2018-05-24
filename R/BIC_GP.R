BIC_GP <- function(Ker_X, Ker_Y, beta){
  n2 <- dim(Ker_X)[1]
  res2 <- mean((Ker_Y - Ker_X %*% beta)^2)
  loglik <- -(n2/2) * (log(2 * pi * res2) + 1)
  df <- length(table(as.vector(beta))) - 1
  bic <- -2 * loglik + log(n2) * df
  return(bic)  
}