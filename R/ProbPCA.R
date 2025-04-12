#' @title Probabilistic PCA
#'
#' @description Function to perform (non-phylogenetic) probabilistic PCA. This function is a modification (fork) of Rdimtools::do.ppca .
#' @param x a matrix with traits in columns and observations in rows.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @returns returns a list of class "phylPPCA". See "Details" for more information.
#' @details
#' This function uses the same algorithm as Rdimtools::do.ppca. However, it returns more details from the estimation and computes AIC and AICc.
#'
#' The function returns a list with the following elements. scores: the scores of the principal components; projection: the eigenvectors; sig: the MLE of the error of the model; mle.W: the MLE of the W matrix: varnames: the names of the traits; loglik: the log-likelihood of the estimate. Function also returns AIC, AICc, and BIC for the model.
#'
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @examples
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ppca <- ProbPCA(x = dt, ret_dim = 2)
#' doBiplot(x = ppca, add_margin = 0.3)
#'
#' @importFrom mclust dmvnorm
#' @importFrom stats cov
#' @export
ProbPCA <- function(x, ret_dim = 2){

  if ((!is.numeric(ret_dim))||(ret_dim<1)||(ret_dim>=ncol(x))||is.infinite(ret_dim)||is.na(ret_dim)){
    stop("*do.pca : 'ndim' is a positive integer in [1,#(covariates)).")
  }

  varnames <- colnames(x)
  N <- nrow(x)
  d <- ncol(x)
  q <- as.integer(ret_dim)
  S <- stats::cov(x)
  eigS <- base::eigen(S, TRUE)
  mlsig2 <- (sum(eigS$values[(q+1):d]))/(d-q)
  mlW <- (eigS$vectors[,1:q])%*%(diag(sqrt(eigS$values[1:q] - mlsig2)))
  M <- (t(mlW)%*%mlW)+(diag(ncol(mlW))*mlsig2)
  ## aux.bicgstab_fork is a fork of Rdimtools:::aux.bicgstab version 1.1.2
  SOL <- aux.bicgstab_fork(M, t(mlW), verbose=FALSE)
  ## aux.adjprojection_fork is a fork of Rdimtools:::aux.adjprojection version 1.1.2
  projection <- aux.adjprojection_fork(t(SOL$x))
  result <- list()
  result$scores <- (x%*%projection)
  result$projection <- projection
  result$sig <- mlsig2
  result$mle.W <- mlW
  result$varnames <- varnames
  ## Likelihood of the model
  muhat <- apply(X = x, MARGIN = 2, FUN = mean)
  S_mle <- mlW%*%t(mlW)+mlsig2*diag(d)
  lik <- sum( mclust::dmvnorm(data = x, mean = muhat, sigma = mlW%*%t(mlW)+mlsig2*diag(d)
                              , log = TRUE) )
  result$loglik <- lik
  params <- (d*q) - (0.5*q*(q-1)) + 1
  result$AIC <- (-2*lik) + (2*params)
  result$AICc <- (-2*lik) + (2*params*(N/(N-params-1)))
  result$BIC <- (-2*lik) + (params*log(N))
  class(result) <- c(class(result), "PPCA")
  return(result)
}
