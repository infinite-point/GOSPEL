\name{GOSPEL}
\alias{GOSPEL}
\title{Graph Oriented Sparse Learning}
\description{Sparse regression for network graphs.}
\usage{
GOSPEL(Ker_X, Ker_Y, CKA_mat, lambda_1 = 0.01, tau, lambda_2)
}
\arguments{
  \item{Ker_X}{Vectorized kernel matrices corresponding to the predictor network graphs.}
  \item{Ker_Y}{A vectorized kernel matrix corresponding to the response network graph.}
  \item{CKA_mat}{A correlation matrix based on the centered kernel alignment.}
  \item{lambda_1}{Candidates of the first regularization parameter in Equation (1).}
  \item{tau}{Candidates of of the threshold for R.}
  \item{lambda_2}{Candidates of the second regularization parameter in Equation (1).}
}
\details{
We propose a sparse learning algorithm for network graph data, called Graph-Oriented SParcE Learning (GOSPEL), 
to find a subset of the topological information matrices of the predictor variables (networks) related to the response variable (network). 
More specifically, we propose to use particular forms of diffusion kernel-based centered kernel alignment as a measure of statistical correlation between graph Laplacian matrices, 
and solve the optimization problem with a novel graph-guided generalized fused lasso.}
\value{
  \item{beta}{A regression coefficient vector.}
  \item{bic}{The score of BIC corresponding to the coefficient vector.}
  \item{lambda_1}{The optimized value of lambda_1.}
  \item{lambda_2}{The optimized value of lambda_2.}
  \item{R}{The reconstructed correlation matrix based on tau.}
  \item{tau}{The optimized value of the threshold.}
}
\author{
Hideko Kawakubo & Teppei Shimamura
}
\references{
Hideko Kawakubo, Yusuke Matsui, Itaru Kushima, Norio Ozaki, and Teppei Shimamura, A Network of Networks Approach for Modeling Interconnected Brain Tissue-Specific Networks, submitted.
}
\examples{
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
}
