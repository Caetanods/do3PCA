
#' Compute the MLE and AIC scores for a model of trait evolution
#'
#' Function will compute the parameter values and the AIC for a model of trait evolution. Models are estimated using a single parameter to model all traits in the multivariate data (x). The models are fit with the product of the likelihood of each independent trait fitted to a tree with the branch lengths transformed under the model M and a single model parameter p. The model parameter p represents an average of the process across all traits. The parameters for the models are: 1) lambda for Pagel's lambda, 2) alpha for OU (with theta as the root value), 3) b for the EB, and 4) sigma for Brownian motion. Please note that fitting BM does not require transformation of branch lengths.
#'
#' @param x a dataset with traits in columns and species on the rows. Please make sure that dataset matches the phylogeny.
#' @param phy a phylogenetic tree with the same species as in x.
#' @param model the name of a model of trait evolution. The options are "BM", "lambda", "OU", or "EB".
#'
#' @returns a list with the elements: model (the model name), lik (the log-likelihood), AIC (the Akaike Information Criterion), log_par (the log of the parameter value), nlopt_out (the output object from nloptr optimizer). If the model is "BM", then the output has R (the evolutionary R matrix).
#' @export
#'
#' @importFrom phytools phyl.vcv rescale
#' @importFrom ratematrix likelihoodFunction
#' @importFrom ape vcv.phylo
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom nloptr nloptr
traitModelAIC <- function(x, phy, model){

  ## Make sure the data has species names:
  if( is.null( rownames(x) ) ){
    rownames(x) <- phy$tip.label
  }

  ## Get minimum value for tolerance of singularity check:
  min_tol <- .Machine$double.xmin * 10

  ## Evolutionary matrix, phy covariance, and phy mean:
  model <- match.arg(arg = model, choices = c("BM", "lambda", "OU", "EB"))

  ## Set some constants:
  N <- nrow(x)
  d <- ncol(x)
  K <- (d * (d + 1)) / 2 ## Number of parameter for the R matrix.

  ## For the BM model we do not need to make branch length transformations:
  if( model == "BM" ){
    ## Just need to get the likelihood here:
    C <- ape::vcv.phylo(phy)
    tt <- try(phytools::phyl.vcv(X = x, C = C, lambda = 1), silent = TRUE)
    lik <- ratematrix::likelihoodFunction(data = x, phy = phy, R = tt$R, root = tt$alpha)
    AIC <- (-2*lik) + 2*K
    ll <- list(model = model, lik = lik, R = tt$R, AIC = AIC)
    return( ll )
  }

  ## Transform branch lengths according to the model:
  if( model == "lambda" ){
    ## Prepare estimation function for lambda
    phy_rescale_fn <- phytools::rescale(x = phy, model = "lambda")
    model_likfn <- function(x){
      l <- x ## Fix for RCRAN.
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
      a <- x ## Fix for RCRAN.
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
      a <- x ## Fix for RCRAN
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

  lik_fn <- function(pars, N = N, d = d){
    pars <- exp( pars ) ## Search in log space.
    branch_trans_eval <- model_likfn(pars)
    res <- branch_trans_eval$lik ## Returning -1 * lik here.
    if( is.infinite(res) ) return(Inf)
    if( is.na(res) ) return(Inf)
    return( res )
  }

  ## Starting states for the parameters:
  pars <- log( runif(n = 1, min = search_bound[1], max = search_bound[2]) )
  ## Check starting point:
  start_p_lik <- lik_fn(pars = pars, N = N, d = d)
  if( is.infinite( start_p_lik ) ) stop("Starting point MLE is infinite!")
  if( is.na( start_p_lik ) ) stop("Starting point MLE is NA!")
  ## The bounds for the search:
  set_lb <- log(search_bound[1])
  set_ub <- log(search_bound[2])
  ## Optimize using nlopt:
  mle <- nloptr::nloptr(x0 = pars, eval_f = lik_fn
                        , lb = set_lb, ub = set_ub
                        , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                      , xtol_rel = 1e-8
                                      , maxeval = 1e6)
                        , N = N, d = d)
  ## Compute the AIC:
  lik <- -1 * mle$objective
  AIC <- (-2*lik) + 2*(K+1) ## Transformation models have +1 parameter.
  ## Return the solution:
  return( list(model = model, lik = lik, AIC = AIC, log_par = mle$solution, nlopt_out = mle) )
}
