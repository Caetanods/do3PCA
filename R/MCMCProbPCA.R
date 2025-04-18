#' @title MCMC for Probabilistic PCA
#'
#' @description Function to perform probabilistic PCA using Bayesian Markov-Chain Monte Carlo.
#' @param x a matrix with traits in columns and species values in rows.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @param gamma_shape The shape parameter for the gamma prior.
#' @param gamma_rate The rate parameter for the gamma prior.
#' @param gen The number generations for the MCMC.
#' @param step_scale The scaler for the mulplier proposal of the MCMC.
#' @param burn The proportion of the burn-in to the excluded before computing the average posterior parameter value.
#' @returns returns a list with "mean_scores", "percent_var_pcs", "var_loadings", and "stats". Also return "mcmc" which is the output from the "mcmc" R package.
#' @details
#' This function uses the "mcmc" R package to conduct a Metropolis-Hastings search. The output "mcmc" can be used for further analyses of the MCMC chain using the "mcmc" package.
#'
#' @references Revell, L. J. 2009. Size-Correction and Principal Components for Interspecific Comparative Studies. Evolution 63:3258–3268. doi: 10.1111/j.1558-5646.2009.00804.x
#' @references Revell, L. J. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ 12:e16505. doi: 10.7717/peerj.16505
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @references Zhang, A., Chan, K.L., Kwok, J.T., and D-Y Yeung. 2004. Proceedings of the AAAI Conference on Artificial Intelligence 19 Book One, Volume. https://aaai.org/papers/00372-aaai04-060-bayesian-inference-on-principal-component-analysis-using-reversible-jump-markov-chain-monte-carlo/
#' @importFrom stats rgamma dgamma
#' @importFrom mcmc metrop
#' @examples
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ## Do a VERY small chain to test:
#' mcmcPPCA <- MCMCProbPCA(x = dt, gen = 100)
#'
#'
#' @export
MCMCProbPCA <- function(x, ret_dim = 2, scale_data = TRUE, gamma_shape = 1, gamma_rate = 2, gen = 10^6, step_scale = 0.08, burn = 0.2){
  ## Standardized data:
  if( scale_data ) x <- scale(x = x)

  ## Preparing parameters:
  d <- ncol(x)
  q <- ret_dim
  N <- nrow(x)
  S <- cov(x) ## Using the MLE of the covariance matrix.
  S_decom <- eigen(S)
  g <- S_decom$values

  ## Starting points:
  L_start <- rgamma(n = q, shape = gamma_shape, rate = gamma_rate)
  sigma_start <- rgamma(n = 1, shape = gamma_shape, rate = gamma_rate)

  ## The MCMC requires the product of the prior and the likelihood:
  post_factory <- function(N, d, g, shape, rate) function(par){
    sigma <- exp(par[length(par)])
    L <- exp(par[1:(length(par)-1)])
    loglikeval <- do3PCA:::loglik(L = L, sigma = sigma, N = N, d = d, g = g)
    prioreval <- sum(dgamma(x = c(L, sigma), shape = gamma_shape
                            , rate = gamma_rate, log = TRUE))
    if( is.na(loglikeval) ){
      print( paste0("NA : ", par ) )
      loglikeval <- Inf
    }
    if( is.infinite(loglikeval) ){
      print( paste0("Inf : ", par ) )
      loglikeval <- Inf
    }
    return( loglikeval + prioreval)
  }
  post_fn <- post_factory(N = N, d = d, g = g, shape = gamma_shape
                          , rate = gamma_rate)

  ## Do the MCMC search:
  print( "Starting MCMC...")
  mcmc1 <- metrop(obj = post_fn, initial = log( c(L_start, sigma_start) )
                  , nbatch = gen, scale = step_scale)
  print( "Done!" )

  burn_in <- floor(gen * burn)
  post <- mcmc1$batch[-1 * (1:burn_in),]

  ## Mean of the posterior for each parameter:
  L_mean_p <- exp( apply(X = post[,1:2], MARGIN = 2, FUN = mean) )
  sigma_mean_p <- exp( mean(post[,3]) )

  ## Reconstruction of W matrix from L vector:
  ## With the W matrix and sigma we can get the principal components.
  L_q <- diag(x = L_mean_p)
  U_q <- S_decom$vectors[,1:q]
  I_q <- diag(nrow = q, ncol = q)
  W <- U_q %*% (L_q - (sigma_mean_p^2*I_q))^0.5 %*% I_q

  ## Compute the projection:
  mcmc_proj <- getProjection(x = x, W = W, S = S, sigma = sigma_mean_p
                             , squared = FALSE)

  ## Compute the percent explained variance of the axes:
  mcmc_pc_var <- mcmc_proj$e_values / sum(mcmc_proj$e_values) * 100

  ## Compute the loadings:
  mcmc_loadings <- mcmc_proj$projection^2 * 100
  rownames(mcmc_loadings) <- mcmc_proj$varnames
  colnames(mcmc_loadings) <- paste0("loadings_PC", 1:ret_dim)

  ## Add to the projection list to return as stats:
  mcmc_proj$mle.W <- NULL ## Change the name to not be confusing!
  mean_scores <- mcmc_proj$scores
  mcmc_proj$scores <- NULL ## Avoid duplicated information!
  mcmc_proj$sigma_mean <- mcmc_proj$sig ## Change name.
  mcmc_proj$sig <- NULL
  mcmc_proj$W_mean <- W
  mcmc_proj$L_mean <- L_mean_p

  ## Add colnames:
  colnames(mean_scores) <- paste0("PC_", 1:ret_dim)

  ## Output of the function:
  out <- list(mean_scores = mean_scores, percent_var_pcs = mcmc_pc_var
              , var_loadings = mcmc_loadings, mcmc = mcmc1, stats = mcmc_proj)
  return( out )
}
