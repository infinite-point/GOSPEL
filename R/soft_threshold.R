soft_threshold <- function(vv, lambda_2){
  n <- nrow(vv)
  p.ind <- which(vv > lambda_2)
  n.ind <- which(vv < - lambda_2)
  res <- Matrix(0, n, 1, sparse=TRUE)
  res[p.ind] <- vv[p.ind] - lambda_2
  res[n.ind] <- vv[n.ind] + lambda_2
  return(res)
}