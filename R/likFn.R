#' @title Probabilistic PCA likelihood function
#'
#' @description Likelihood function for the standard probabilistic PCA or the phylogenetic probabilistic PCA. Note that one just need to change between the covariance of the data and the R matrix (multivariate Brownian Motion covariance matrix). The likelihood function only uses the eigenvalues of the covariance matrix (see parameter g).
#' @param L vector with length equal to the number of retained dimensions.
#' @param sigma variance. This is the average variance lost per discarded dimension.
#' @param N sample size. Number of observations.
#' @param d integer. Number of dimensions in the data.
#' @param g eigenvalues of the covariance matrix.
#' @returns returns the -1 * log-likelihood of the data.
#' @references Revell, L. J. 2009. Size-Correction and Principal Components for Interspecific Comparative Studies. Evolution 63:3258–3268. doi: 10.1111/j.1558-5646.2009.00804.x
#' @references Revell, L. J. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ 12:e16505. doi: 10.7717/peerj.16505
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @references Zhang, A., Chan, K.L., Kwok, J.T., and D-Y Yeung. 2004. Proceedings of the AAAI Conference on Artificial Intelligence 19 Book One, Volume. https://aaai.org/papers/00372-aaai04-060-bayesian-inference-on-principal-component-analysis-using-reversible-jump-markov-chain-monte-carlo/
#' @examples
#' \donttest{
#' phy <- ratematrix::anoles$phy[[1]]
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' L_start <- rgamma(n = 2, shape = 2, rate = 1)
#' sigma_start <- rgamma(n = 1, shape = 2, rate = 1)
#' R <- do3PCA:::fastR(x = dt, phy = phy)
#' g <- eigen(R)$values
#' loglik(L = L_start, sigma = sigma_start, N = nrow(dt), d = 2, g = g)
#' }
#'
#' @export
loglik <- function(L, sigma, N, d, g){
  q <- length(L) # Dimension of the model.
  pA <- ((N*d)/2) * log(2*pi)
  pB <- sum( sapply(1:q, function(i) (N/2)*log(L[i]) + N*(d-q)*log(sigma)) )
  pC <- N/2 * sum( sapply(1:q, function(i) (L[i]^(-1)) * g[i]) )
  pD <- (N * sigma^(-2))/2 * sum(g[(q+1):d])
  return( -1 * (pA+pB+pC+pD) )
}
