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