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