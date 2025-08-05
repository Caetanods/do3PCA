# LaplaceApprox <- function(lambda, N){
#   ## Function to compute the automatic number of principal components as described in Minka. https://tminka.github.io/papers/pca/
#   d <- length(lambda) ## Total number of dimensions/features
#   k_vec <- 2:(d-1) ## Models considered
#   pD <- vector(mode = "numeric", length = length(k_vec)) ## p(D|k)
#   names(pD) <- paste0(k_vec, " PCs")
#
#   ## ###########################################################################
#   ## Block to compute log(|Az|). Part of code from: https://github.com/WeiAkaneDeng/SPAC2 written by Wei Deng under GPL-2 license. Modified on July 23, 2025.
#   logediff <- vector(mode = "numeric", length = length(lambda)-1)
#   for (i in 1:(length(lambda)-1)){
#     j <- (i+1):length(lambda)
#     logediff[i] <- sum(log(lambda[i] - lambda[j]))
#   }
#   cumsum_logediff <- cumsum(logediff)
#   inve <- 1/lambda
#   invediff <- matrix(rep(inve, length(lambda)), ncol=length(lambda)) - matrix(rep(inve, length(lambda)),nrow=length(lambda), byrow=T)
#   invediff[invediff <= 0 ] <- 1 ## The entire upper-diagonal is negative!
#   invediff <- log(invediff) ## Upper diagonal is all zeros!
#   ## Note that next we have a sum, so 0s might be on purpose.
#   cumsum_invediff <- apply(invediff, 2, cumsum)
#   row_invediff <- t(apply(cumsum_invediff, 1, cumsum))
#   ## END Block
#   ## ###########################################################################
#
#   ## Make the Laplace Approximation for each model
#   for( k in k_vec ){
#     v <- mean(lambda[(k+1):d]) # variance/error
#     ## Compute log(|Az|), see block above.
#     lAz <- row_invediff[k] + cumsum_logediff[k]
#     lAz <- lAz + (d-k) * sum(log(1/v - inve[1:k]))
#     m <- (d*k) - (k*(k+1)/2)
#     lAz <- lAz + m*log(N)
#     ## Compute log(p(U))
#     sumGamma <- sapply(1:k, function(x) lgamma((d-x+1)/2)-((d-x+1)/2*log(pi)))
#     sumGamma <- sum(sumGamma)
#     lpU <- -k*log(2) + sumGamma
#     ## Get the terms for the final sum (see Minka equation 30, log-transformed)
#     t2 <- (-N/2) * sum( log(lambda[1:k]) )
#     t3 <- ((-N*(d-k))/2) * log(v)
#     t4 <- ((m+k)/2) * log(2*pi)
#     t5 <- -0.5 * lAz
#     t6 <- -k/2 * log(N)
#     pD[k-1] <- lpU + t2 + t3 + t4 + t5 + t6 ## Minka equation 30 in log space
#   }
#   ## Best model:
#   nPC <- unname( which.max(pD)+1 ) # Note that starts from 2
#   return( list("best # PCs" = nPC, "log(p(D|k))" = pD) )
# }

#' @noRd
#' @importFrom stats setNames
LaplaceApprox <- function(lambda, N, tau=0.001) {
  ## Function to compute the automatic number of principal components as described in Minka. https://tminka.github.io/papers/pca/
  ## Code adapted from: https://github.com/WeiAkaneDeng/SPAC2 written by Wei Deng under GPL-2 license. Modified on August 5, 2025.

  ## Drop eigenvalues that are very small. Small than the min threshold (tau)
  lambda <- lambda[lambda > tau]
  n <- length(lambda) ## Number of actual dimensions, after dropping small eigenvalues (if present).
  kmax <- min(c(n-1, N-2))

  ## Computing the difference in the log.
  logediff <- NA
  for (i in 1:(length(lambda)-1)){
    j <- (i+1):length(lambda)
    ## Second term in the sum will only exist if there is a very small eigenvalue in the lambda vector.
    logediff[i] <- sum(log(lambda[i] - lambda[j])) + (n-length(lambda))*log(lambda[i])
  }
  ## The cumulative sum of the log difference between the eigenvalue_i and all later eigenvalues.
  cumsum_logediff <- cumsum(logediff)
  n1 <- N-1 ## N - 1, similar term used when computing the sample variance.
  ehat <- (lambda)*N/n1
  inve <- 1/lambda
  ## Each row is the original inve vector minus the ith element of the original inve vector. A quick way to perform the subtraction. This is avoiding a loop.
  invediff <- matrix(rep(inve, length(lambda)), ncol=length(lambda)) - matrix(rep(inve, length(lambda)),nrow=length(lambda), byrow=T)
  invediff[invediff <= 0]<- 1 ## The entire upper-diagonal is negative!
  invediff <- log(invediff) ## Upper diagonal is all zeros!
  ## This is the same as the previous loop. This is the same as "logediff" above. Using a matrix trick which is neat.
  cumsum_invediff <- apply(invediff, 2, cumsum)
  row_invediff <- t(apply(cumsum_invediff, 1, cumsum))

  loge <- log(ehat)
  cumsum_loge <- cumsum(loge)
  cumsum_e <- cumsum(ehat)

  dn <- length(lambda)
  kmax <- length(lambda)-1
  ks <- 1:kmax
  z <- (n-ks+1)/2*log(pi) - lgamma((n-ks+1)/2)
  ## Given the notes below, cumsum_z is -log(p(U))
  cumsum_z <- cumsum(z)
  ## Here is the fix added:
  cumsum_z <- cumsum_z + log(2)*ks

  ## This is the computation of the probability of the data given the model.
  ## This is where we get the value that we use for the model selection.
  ## This should match equation 76 and 77 of the Minka Report OR equations 30 and 31 of the paper.

  prb <- NA
  for (i in 1:length(ks)){
    k <- i
    # v is computed as the sum of all lambda minus the lambda in the model. A trick for the computation.
    v <- (cumsum_e[length(cumsum_e)] - cumsum_e[k])/(n-k) # the variance
    ## cumsum_loge is the sum of the eigenvalues in part 2
    prb[i] <- -n1/2*cumsum_loge[k] + (-n1*(n-k)/2)*log(v) # part 2 + part 3
    prb[i] <- prb[i] - cumsum_z[k] - k/2*log(n1) # part 1 + part 6
    h <- row_invediff[k] + cumsum_logediff[k]
    h <- h + (n-k)*sum(log(1/v - inve[1:k]))
    mm <- n*k-k*(k+1)/2
    h <- h + mm*log(N)
    prb[i] <- prb[i] + (mm+k)/2*log(2*pi) - h/2 # part 4 + part 5
    prb[i] <- prb[i] + 1.5*k*log(2)
    prb[i] <- prb[i] - 0.5*log(n-k)
  }
  prb <- setNames(object = prb, nm = 1:length(prb))
  return( list("best # PCs" = which.max(prb), "log(p(D|k))" = prb) )
}

#' Automatic choice of number of principal components
#'
#' @param x data. Columns are variables/features. Rows are observations/species.
#' @param phy phylogenetic tree. Number of tips must match number of rows in x.
#' @param cov a covariance matrix. Either the covariance of the data or an evolutionary ratematrix (R).
#' @param N number of observations/species in the tree.
#'
#' @returns a list with the best fitting number of principal components (k) and a vector with the log probability of the data given the number of principal components (k).
#' @export
#' @importFrom ape vcv.phylo
modelChoice <- function(x = NULL, phy = NULL, cov = NULL, N = NULL){
  ## Function works with these combinations of parameters:
  ## x; model choice for standard pPCA with covariance estimation
  ## x and phy; model choice for 3PCA with covariance estimation
  ## cov and N; model choice based on provided covariance (standard or phylogenetic) and the number of samples or species in the phylogeny.

  has_x <- !is.null(x)
  has_phy <- !is.null(phy)
  has_cov <- !is.null(cov)
  has_N <- !is.null(N)

  if(has_x & has_phy & has_cov & has_N){
    message("All arguments present! Function just needs x OR x and phy OR cov and N. Please check help page.")
  }

  if(has_x & !has_phy & !has_cov){
    ## Standard pPCA
    lambda <- eigen( cov(x) )$val
    N <- nrow(x)
  } else if(has_x & has_phy & !has_cov){
    ## Phylogenetic pPCA (Here assuming BM model. Can extend this later.)
    C <- ape::vcv.phylo(phy)
    invC <- solve(C)
    a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
    A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
    R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
    lambda <- eigen( R )$val
    N <- nrow(x)
  } else if(has_cov & has_N){
    ## Use the cov matrix
    lambda <- eigen( cov )$val
  } else{
    stop("Incorrect inputs. Please check help page.")
  }

  return( LaplaceApprox(lambda = lambda, N = N) )
}
