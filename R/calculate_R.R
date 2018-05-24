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