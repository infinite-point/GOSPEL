# GOSPEL
R package for graph-oriented sparce learning

<strong>Depends:</strong>

R(>=3.4.1)

data.table, rARPACK, igraph, bayesopt, Matrix

<strong>Authors:</strong>

Hideko Kawakubo and Teppei Shimamura

<strong>Contact:</strong>
kawakubo[at]med.nagoya-u.ac.jp and shimamura[at]med.nagoya-u.ac.jp

## General overview

We propose a sparse learning algorithm for network graph data, called Graph-Oriented SParcE Learning (GOSPEL), to find a subset of the topological information matrices of the predictor variables (networks) related to the response variable (network). More specifically, we propose to use particular forms of diffusion kernel-based centered kernel alignment as a measure of statistical correlation between graph Laplacian matrices, and solve the optimization problem with a novel graph-guided generalized fused lasso.

## An example

```
library(GOSPEL)
data(additive_random_500_XY)
p <- dim(X)[3]
gamma_list <- 10^seq(-1, 1, length=5)
CKA_mat <- array(0, dim=c(p,p,length(gamma_list)))
result_mat <- matrix(0, p, length(gamma_list))
BIC_mat <- c(0, length(gamma_list))
for (i in 1:length(gamma_list)){
  cat("gamma list:",i,"\n")
  gamma <- gamma_list[i]
  obj <- calculate_R(X, Y, gamma)
  CKA_mat[,,i] <- obj$CKA_mat
  Ker_X <- obj$Ker_X
  Ker_Y <- obj$Ker_Y
  ans <- GOSPEL(Ker_X, Ker_Y, CKA_mat[,,i])
  result_mat[,i] <- as.matrix(ans$beta)
  BIC_mat[i] <- ans$bic
}
min_id <- which.min(BIC_mat)
beta <- result_mat[,min_id]
beta
```

## Reference
Hideko Kawakubo, Yusuke Matsui, Itaru Kushima, Norio Ozaki, and Teppei Shimamura, A Network of Networks Approach]{A Network of Networks Approach for Modeling Interconnected Brain Tissue-Specific Networks, submitted.
