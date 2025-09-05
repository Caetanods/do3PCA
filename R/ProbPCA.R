#' @title Probabilistic PCA
#'
#' @description Function to perform (non-phylogenetic) probabilistic PCA.
#' @param x a matrix with traits in columns and observations in rows.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @param scale if data should be rescaled to have zero mean and unit variance using the function scale(). Default is TRUE.
#' @returns returns a list.
#' @details
#' Optimization is performed using nloptr.
#'
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @examples
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ppca <- ProbPCA(x = dt, ret_dim = 2)
#'
#' @importFrom mclust dmvnorm
#' @importFrom stats cov
#' @importFrom nloptr nloptr
#' @export
ProbPCA <- function(x, ret_dim = 2, scale = TRUE){

  if( !inherits(x = x, what = "matrix") ){
    stop("x needs to be a matrix. Cannot be a data.frame.")
  }

  if( scale ){
    ## Perform tha scale:
    x <- scale(x)
  }

  if ((!is.numeric(ret_dim))||(ret_dim<2)||(ret_dim>=ncol(x))||is.infinite(ret_dim)||is.na(ret_dim)){
    d <- ncol(x)
    stop( paste0("'ret_dim' must be a positive integer between 2 and ", d-1, ".", collapse = "") )
  }
  ## Compute quantities:
  varnames <- colnames(x)
  N <- nrow(x)
  d <- ncol(x)
  q <- as.integer(ret_dim)
  S <- stats::cov(x)
  S_decom <- eigen(S)
  g <- S_decom$values
  ## Likelihood function in log parameter space:
  log_space <- function(pars, N, d, g){
    ll <- loglik(L = exp(pars[1:q]), sigma = exp(pars[q+1])
                 , N = N, d = d, g = g)
    if( is.infinite(ll) ) ll <- NA
    if( is.na(ll) ){
      return( Inf ) ## NLOPT is minimizing!
    } else{
      return( -1 * ll ) ## NLOPT is minimizing!
    }
  }
  ## Starting point for the optimization:
  pars <- log( rgamma(n = q+1, shape = 1, rate = 2) )
  ## Check starting point:
  start_p_lik <- log_space(pars = pars, N = N, d = d, g = g)
  if( is.infinite( start_p_lik ) ) stop("Starting point MLE is infinite!")
  if( is.na( start_p_lik ) ) stop("Starting point MLE is NA!")
  ## Optimize using nlopt:
  mle <- nloptr::nloptr(x0 = pars, eval_f = log_space
                        , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                      , xtol_rel = 1e-8
                                      , maxeval = 1e6)
                        , N = N, d = d, g = g)
  ## MLE solution:
  L_mle <- exp(mle$solution[1:q])
  sigma_mle <- exp(mle$solution[q+1])
  best_lik <- -1 * mle$objective
  ## Reconstruction of W matrix from L vector:
  ## With the W matrix and sigma we can get the principal components.
  L_q <- diag(x = L_mle)
  U_q <- S_decom$vectors[,1:q]
  I_q <- diag(nrow = q, ncol = q)
  W <- U_q %*% (L_q - (sigma_mle^2*I_q))^0.5 %*% I_q
  ## Compute the projection:
  mle_proj <- getProjection(x = x, W = W, S = S, sigma = sigma_mle
                            , squared = FALSE)
  ## Compute the percent explained variance of the axes:
  mle_pc_var <- mle_proj$e_values / sum(mle_proj$e_values) * 100
  ## Compute the loadings:
  mle_loadings <- mle_proj$projection^2 * 100
  rownames(mle_loadings) <- mle_proj$varnames
  colnames(mle_loadings) <- paste0("loadings_PC", 1:ret_dim)

  ## Add to the projection list to return as stats:
  scores <- mle_proj$scores
  mle_proj$scores <- NULL ## Avoid duplicated information!
  mle_proj$mle.L <- L_mle
  mle_proj$loglik <- best_lik
  params <- (d*q) - (0.5*q*(q-1)) + 1
  mle_proj$npar <- params

  ## Not using AIC and BIC.
  # mle_proj$AIC <- (-2*best_lik) + (2*params)
  # mle_proj$AICc <- (-2*best_lik) + (2*params*(N/(N-params-1)))
  # mle_proj$BIC <- (-2*best_lik) + (params*log(N))

  ## Define the output and class:
  out_obj <- list(mle_scores = scores, percent_var_pcs = mle_pc_var
                  , var_loadings = mle_loadings, log_lik = mle_proj$loglik
                  , nlopt = mle, stats = mle_proj)
  class( out_obj ) <- append(x = class( out_obj), values = "PPCA")

  ## Need to create class for the output object!
  return( out_obj )
}
