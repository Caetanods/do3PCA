#' @title Probabilistic Phylogenetic PCA (Joint)
#'
#' @description Function to perform probabilistic phylogenetic PCA. Allows for fit of alternative models of trait evolution using branch length transformation.
#' @param phy phylogeny in "phylo" format.
#' @param x a matrix with traits in columns and species values in rows. Rownames must match the tip labels of phylogeny.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @param model choice of model of trait evolution. One of "BM", "lambda", "OU", or "EB".
#' @param scale if data should be rescaled to have zero mean and unit variance using the function scale(). Default is TRUE.
#' @param quiet if function should suppress output to the console while running
#' @returns returns a list of class "phylPPCA". See "Details" for more information.
#' @details
#' The function can be used to estimate the probabilistic phylogenetic PCA (3PCA) using distinct models of trait evolution. Models are implemented using branch length transformation. Model fitting happens in two steps. First the maximum likelihood of the evolutionary covariance matrix (R) and the parameter of the model is estimated. Then the 3PCA model is estimated using the phylogenetic tree with branch lengths transformed following the MLE for the parameter of each trait evolution model.
#'
#' The function returns a list with the following elements. scores: the scores of the principal components; e_values: eigenvalues; e_vectors: eigenvectors or the projection; model_fit: information about the trait evolution model; loadings: the loadings of the PCs; varnames: the names of the variables; sig: the MLE of the error; mle.W: the MLE of the W matrix; Function also returns AIC, AICc, and BIC for the model.
#' @references Revell, L. J. 2009. Size-Correction and Principal Components for Interspecific Comparative Studies. Evolution 63:3258–3268. doi: 10.1111/j.1558-5646.2009.00804.x
#' @references Revell, L. J. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ 12:e16505. doi: 10.7717/peerj.16505
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @importFrom ape vcv.phylo
#' @importFrom phytools rescale phyl.vcv
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom ratematrix likelihoodFunction
#' @importFrom nloptr nloptr
#' @importFrom stats lm coef
#' @importFrom stats setNames runif
#' @importFrom Matrix solve
#' @examples
#' phy <- ratematrix::anoles$phy[[1]]
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ppca <- phylProbPCA(phy = phy, x = dt, ret_dim = 2)
#' doBiplot(x = ppca, add_margin = 0.3)
#'
#' @export
phylProbPCAJoint <- function(phy, x, ret_dim = 2, model = "BM", scale = TRUE, quiet = FALSE){
  ## Check if x is a matrix object.
  if( !inherits(x = x, what = "matrix") ){
    stop("x needs to be a matrix. Cannot be a data.frame.")
  }

  if( scale ){
    ## Perform tha scale:
    x <- scale(x)
  }

  ## Check parameters
  if ((!is.numeric(ret_dim))||(ret_dim<2)||(ret_dim>=ncol(x))||is.infinite(ret_dim)||is.na(ret_dim)){
    stop("ret_dim is an integer between 2 and the number of variables -1.")
  }

  ## Get minimum value for tolerance of singularity check:
  min_tol <- .Machine$double.xmin * 10

  ## Evolutionary matrix, phy covariance, and phy mean:
  model <- match.arg(arg = model, choices = c("BM", "lambda", "OU", "EB"))

  ## Set some constants:
  N <- nrow(x)
  d <- ncol(x)
  q <- as.integer(ret_dim)

  ## Transform branch lengths according to the model:
  if( model == "BM" ){
    ## Branch lengths remain the same so likelihood is fixed:
    C <- ape::vcv.phylo(phy)
    ## Function can break if C cannot be inverted:
    tt <- try(phytools::phyl.vcv(X = x, C = C, lambda = 1), silent = TRUE)
    if( inherits(x = tt, what = "try-error") ){
      lik_BM <- Inf ## Negative log-lik
    } else{
      g_BM <- eigen(tt$R)$values
      lik_BM <- -1 * ratematrix::likelihoodFunction(data = x, phy = phy, R = tt$R, root = tt$alpha)
    }
  }
  if( model == "lambda" ){
    ## Prepare estimation function for lambda
    phy_rescale_fn <- phytools::rescale(x = phy, model = "lambda")
    model_likfn <- function(x){
      l <- x # Correct CRAN Note.
      phy_lambda <- phy_rescale_fn(lambda = l)
      ## Matrix can be invertible even if branch lengths are negative
      if( any( phy_lambda$edge.length <= 0 ) ){
        return( list(lik = Inf, R = NA) )
      }
      C <- ape::vcv.phylo(phy_lambda)
      tt <- try(phytools::phyl.vcv(X = x, C = C, lambda = 1), silent = TRUE)
      if( inherits(x = tt, what = "try-error") ){
        return( list(lik = Inf, R = NA) )
      } else{
      ## Inverse of the loglik (nloptr is minimizing!)
      return( list(lik = -1 * ratematrix::likelihoodFunction(data = x, phy = phy_lambda, R = tt$R
                                                             , root = tt$alpha), R = tt$R) )
      }
    }
    ## Search the bounds for the lambda parameter. This function needs to hide warnings, because search will hit out of bounds.
    fn_to_search_l <- function(x){
      ## TRUE if it is bad.
      phy_lambda <- suppressWarnings( phy_rescale_fn(lambda = x) )
      if( any( phy_lambda$edge.length <= 0 ) ) return(TRUE)
      ## Matrix can be invertible even if branch lengths are negative.
      return( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_lambda), tol = min_tol) )
    }
    bound <- binary_search(fn = fn_to_search_l, lb = 0, ub = 1000, chunk = 0.000001)
    search_bound <- c(0, bound[1])
  }
  if( model == "OU" ){
    ## Prepare estimation function for OU
    ## Need to check this estimation. The goal is to simulate a dataset using a multivariate OU model with alpha as a diagonal matrix but having correlation in the R matrix. Then see if we can estimate alpha first using rescale.
    phy_rescale_fn <- phytools::rescale(x = phy, model = "OU")
    model_likfn <- function(x){
      a <- x # Correct CRAN Note.
      phy_ou <- phy_rescale_fn(alpha = a)
      C <- ape::vcv.phylo(phy_ou)
      tt <- try(phytools::phyl.vcv(X = x, C = C, lambda = 1), silent = TRUE)
      if( inherits(x = tt, what = "try-error") ){
        return( list(lik = Inf, R = NA) )
      } else{
        return( list(lik = -1 * ratematrix::likelihoodFunction(data = x, phy = phy_ou, R = tt$R
                                                               , root = tt$alpha), R = tt$R) )
      }
    }
    ## Find maximum value for alpha:
    fn_to_search_ou <- function(x){
      ## TRUE if it is bad.
      phy_ou <- phy_rescale_fn(alpha = x)
      return( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_ou), tol = min_tol) )
    }
    bound <- binary_search(fn = fn_to_search_ou, lb = 0, ub = 1000, chunk = 0.000001)
    search_bound <- c(0, bound[1])
  }
  if( model == "EB" ){
    ## Note that for EB model the max parameter value is 0. So the parameter space is negative!
    phy_rescale_fn <- phytools::rescale(x = phy, model = "EB")
    model_likfn <- function(x){
      a <- x # Correct CRAN Note.
      phy_eb <- phy_rescale_fn(a = -1 * a) ## Searching negative space.
      C <- ape::vcv.phylo(phy_eb)
      ## Need to have another check for matrix singularity:
      if( matrixcalc::is.singular.matrix(x = C, tol = min_tol) ){
        return( list(lik = Inf, R = NA) )
        }
      tt <- try(phytools::phyl.vcv(X = x, C = C, lambda = 1), silent = TRUE)
      if( inherits(x = tt, what = "try-error") ){
        return( list(lik = Inf, R = NA) )
      } else{
        return( list(lik = -1 * ratematrix::likelihoodFunction(data = x, phy = phy_eb, R = tt$R
                                                               , root = tt$alpha), R = tt$R) )
      }
    }
    ## Search for the minimum value for EB parameter.
    fn_to_search_eb <- function(x){
      phy_eb <- phy_rescale_fn(a = -1 * x) ## Searching negative space.
      ## Need to have another check for matrix singularity:
      return( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_eb), tol = min_tol ) )
    }
    bound <- binary_search(fn = fn_to_search_eb, lb = 0, ub = 1000, chunk = 0.000001)
    search_bound <- c(0, bound[1])
  }

  ## ###########################################################################
  ## Joint likelihood function:
  if( model == "BM" ){
    joint_lik_fn <- function(pars, N = N, d = d){
      pars <- exp( pars ) ## Search in log space.
      ppca_lik <- -1 * loglik(L = pars[1:q], sigma = pars[q+1], N = N, d = d, g = g_BM)
      ## Both likelihoods return -1*lik because NLOPTR is minimizing!
      res <- lik_BM + ppca_lik
      ## Protects both positive and negative Inf
      if( is.infinite(res) ) return(Inf)
      if( is.na(res) ) return(Inf)
      return( res )
    }
  } else{
    joint_lik_fn <- function(pars, N = N, d = d){
      pars <- exp( pars ) ## Search in log space.
      branch_trans_eval <- model_likfn(pars[q+2])
      model_lik <- branch_trans_eval$lik ## Returning -1 * lik here.
      if( is.infinite(model_lik) ) return(Inf)
      if( is.na(model_lik) ) return(Inf)
      g <- eigen(branch_trans_eval$R)$values
      ppca_lik <- -1 * loglik(L = pars[1:q], sigma = pars[q+1], N = N, d = d, g = g)
      ## Both likelihoods return -1*lik because NLOPTR is minimizing!
      res <- model_lik + ppca_lik
      if( is.infinite(res) ) return(Inf)
      if( is.na(res) ) return(Inf)
      return( res )
    }
  }

  ## ###########################################################################

  ## Starting states for the parameters:
  pars <- log( rgamma(n = q+1, shape = 1, rate = 2) )
  ## If the model is not BM, then we need another parameter
  if(model != "BM") pars[q+2] <- log( runif(n = 1, min = search_bound[1], max = search_bound[2]) )
  ## Check starting point:
  start_p_lik <- joint_lik_fn(pars = pars, N = N, d = d)
  if( is.infinite( start_p_lik ) ) stop("Starting point MLE is infinite!")
  if( is.na( start_p_lik ) ) stop("Starting point MLE is NA!")
  ## Set the bounds for the search:
  if(model == "BM"){
    set_lb <- rep(x = -Inf, times = q+1)
    set_ub <- rep(x = Inf, times = q+1)
  } else{
    set_lb <- c(rep(x = -Inf, times = q+1), log(search_bound[1]))
    set_ub <- c(rep(x = Inf, times = q+1), log(search_bound[2]))
  }
  ## Optimize using nlopt:
  mle <- nloptr::nloptr(x0 = pars, eval_f = joint_lik_fn
                        , lb = set_lb, ub = set_ub
                        , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                      , xtol_rel = 1e-8
                                      , maxeval = 1e6)
                        , N = N, d = d)
  ## Compute quantities given the MLE solution:
  if( model == "BM" ){
    phy_scaled <- phy
  } else if( model == "EB"){
    phy_scaled <- phy_rescale_fn(-1 * exp(mle$solution[q+2]))
  } else{
    phy_scaled <- phy_rescale_fn(exp(mle$solution[q+2]))
  }
  C <- ape::vcv.phylo(phy_scaled)
  ## Try to add error protection here:
  p_vcv <- phytools::phyl.vcv(X = x, C = C, lambda = 1)
  #invC <- solve(C)
  invC <- Matrix::solve(C, sparse = TRUE)
  #a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
  #A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
  #R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
  R_decom <- eigen(p_vcv$R)
  ## MLE solution:
  L_mle <- exp(mle$solution[1:q])
  sigma_mle <- exp(mle$solution[q+1])
  if( model == "BM" ){
    model_par_mle <- NA
  } else if(model == "EB"){
    model_par_mle <- -1 * exp(mle$solution[q+2])
  } else{
    model_par_mle <- exp(mle$solution[q+2])
  }
  best_lik <- -1 * mle$objective
  ## Reconstruction of W matrix from L vector:
  ## With the W matrix and sigma we can get the principal components.
  L_q <- diag(x = L_mle)
  U_q <- R_decom$vectors[,1:q]
  I_q <- diag(nrow = q, ncol = q)
  W <- U_q %*% (L_q - (sigma_mle^2*I_q))^0.5 %*% I_q
  ## Compute the projection:
  ## Check if a should be a matrix or vector!
  mle_proj <- getPhyloProjection(x = x, invC = invC, R = p_vcv$R, R_decom = R_decom
                                 , W = W, a = p_vcv$alpha, N = N, d = d, q = q
                                 , sigma = sigma_mle)
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

  ## Number of parameters depend on model:
  if( model == "BM" ){
    params <- (d*q) - (0.5*q*(q-1)) + 1
  } else{
    ## One extra parameter for all other models.
    params <- (d*q) - (0.5*q*(q-1)) + 2
  }
  mle_proj$npar <- params
  mle_proj$AIC <- (-2*best_lik) + (2*params)
  mle_proj$AICc <- (-2*best_lik) + (2*params*(N/(N-params-1)))
  mle_proj$BIC <- (-2*best_lik) + (params*log(N))
  mle_proj$phy_scale_par <- exp(mle$solution[q+2])
  mle_proj$phy_scale_model <- model

  ## ###########################################################################
  ## Block to compute the log-lik for the branch length transformation and the pPCA model in separate. Return these values.
  if( model == "BM" ){
    pPCA_lik <- loglik(L = exp(mle$solution[1:q]), sigma = exp(mle$solution[q+1])
                       , N = N, d = d, g = g_BM)
    model_lik <- -1 * lik_BM
  } else{
    branch_trans_eval <- model_likfn(exp(mle$solution[q+2]))
    model_lik <- -1 * branch_trans_eval$lik ## Returning -1 * lik here.
    g <- eigen(branch_trans_eval$R)$values
    pPCA_lik <- loglik(L = exp(mle$solution[1:q]), sigma = exp(mle$solution[q+1])
                       , N = N, d = d, g = g)
  }
  sep_lik <- setNames(object = c(pPCA_lik, model_lik), nm = c("pPCA", "model"))

  ## Need to create class for the output object!
  return( list(mle_scores = scores, percent_var_pcs = mle_pc_var, var_loadings = mle_loadings
               , log_lik = mle_proj$loglik, sep_log_lik = sep_lik, AICc = mle_proj$AICc
               , nlopt = mle, stats = mle_proj) )
}
