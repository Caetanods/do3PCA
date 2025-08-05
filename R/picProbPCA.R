#' @title Probabilistic PCA of independent contrasts (PIC-PCA)
#'
#' @description Function to perform probabilistic PCA of independent contrasts. Allows for fit of alternative models of trait evolution using branch length transformation.
#' @param phy phylogeny in "phylo" format.
#' @param x a matrix with traits in columns and species values in rows. Rownames must match the tip labels of phylogeny.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @param model choice of model of trait evolution. One of "BM", "lambda", "OU", or "EB".
#' @param scale if data should be rescaled to have zero mean and unit variance using the function scale(). Default is TRUE.
#' @param quiet if function should suppress output to the console while running
#' @returns returns a list of class "phylPPCA". See "Details" for more information.
#' @details
#' The function can be used to estimate PIC-PCA using distinct models of trait evolution. Models are implemented using branch length transformation. Model fitting happens in two steps. First the maximum likelihood of the evolutionary covariance matrix (R) and the parameter of the model is estimated. Then the PIC-PCA model is estimated using the phylogenetic tree with branch lengths transformed following the MLE for the parameter of each trait evolution model.
#'
#' The function returns a list with the following elements. scores: the scores of the principal components; e_values: eigenvalues; e_vectors: eigenvectors or the projection; model_fit: information about the trait evolution model; loadings: the loadings of the PCs; varnames: the names of the variables; sig: the MLE of the error; mle.W: the MLE of the W matrix; Function also returns AIC, AICc, and BIC for the model.
#' @references Revell, L. J. 2009. Size-Correction and Principal Components for Interspecific Comparative Studies. Evolution 63:3258–3268. doi: 10.1111/j.1558-5646.2009.00804.x
#' @references Revell, L. J. 2024. phytools 2.0: an updated R ecosystem for phylogenetic comparative methods (and other things). PeerJ 12:e16505. doi: 10.7717/peerj.16505
#' @references Garland Jr, T., Harvey, P.H. and Ives, A.R., 1992. Procedures for the analysis of comparative data using phylogenetically independent contrasts. Systematic biology, 41(1), pp.18-32.
#' @references Tipping, M. E., and C. M. Bishop. 1999. Probabilistic Principal Component Analysis. Journal of the Royal Statistical Society Series B: Statistical Methodology 61(3):611–622. doi: 10.1111/1467-9868.00196
#' @importFrom ape vcv.phylo pic
#' @importFrom phytools rescale phyl.vcv
#' @importFrom matrixcalc is.singular.matrix
#' @importFrom ratematrix likelihoodFunction
#' @importFrom nloptr nloptr
#' @importFrom stats lm coef
#' @examples
#' phy <- ratematrix::anoles$phy[[1]]
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ppca <- picProbPCA(phy = phy, x = dt, ret_dim = 2)
#' doBiplot(x = ppca, add_margin = 0.3)
#'
#' @export
picProbPCA <- function(phy, x, ret_dim = 2, model = "BM", scale = TRUE, quiet = FALSE){

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
    ## The standard implementation.
    pic <- apply(X = x, MARGIN = 2, FUN = pic, phy = phy)
  }
  if( model == "lambda" ){
    ## Prepare estimation function for lambda
    phy_lambda_fn <- phytools::rescale(x = phy, model = "lambda")
    likfn_l <- function(l){
      phy_lambda <- phy_lambda_fn(lambda = l)
      ## Matrix can be invertible even if branch lengths are negative
      if( any( phy_lambda$edge.length <= 0 ) ) return( Inf )
      C <- ape::vcv.phylo(phy_lambda)
      ## Check if the matrix is singular
      if( matrixcalc::is.singular.matrix(x = C, tol = min_tol) ) return( Inf )
      tt <- phytools::phyl.vcv(X = x, C = C, lambda = 1)
      ## Inverse of the loglik (nloptr is minimizing!)
      return( -1 * ratematrix::likelihoodFunction(data = x, phy = phy_lambda, R = tt$R, root = tt$alpha) )
    }
    ## Search the bounds for the lambda parameter. This function needs to hide warnings, because search will hit out of bounds.
    fn_to_search_l <- function(x){
      ## TRUE if 0it is bad.
      phy_lambda <- suppressWarnings( phy_lambda_fn(lambda = x) )
      if( any( phy_lambda$edge.length <= 0 ) ) return(TRUE)
      ## Matrix can be invertible even if branch lengths are negative.
      if( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_lambda), tol = min_tol) ){
        return(TRUE)
      } else{
        return(FALSE)
      }
    }
    bound <- binary_search(fn = fn_to_search_l, lb = 0, ub = 1000, chunk = 0.000001)
    ## Estimate MLE for lambda parameter:
    fit_l <- nloptr::nloptr(x0 = bound[1]/2, eval_f = likfn_l, lb = 0.0, ub = bound[1]
                            , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                          , xtol_rel = 1e-04, print_level = 0))
    ## This is a vector with the model specifics
    fit_model_solution <- list(model = "lambda", lambda = fit_l$solution, npar = 1
                               , loglik = -1*fit_l$objective)
    ## Use the MLE of the lambda to compute the probabilistic PCA model.
    phy_lambda <- phy_lambda_fn(lambda = fit_model_solution[[2]])
    ## Compute the pics with the rescaled tree:
    pic <- apply(X = x, MARGIN = 2, FUN = pic, phy = phy_lambda)
  }
  if( model == "OU" ){
    ## Prepare estimation function for OU
    ## Need to check this estimation. The goal is to simulate a dataset using a multivariate OU model with alpha as a diagonal matrix but having correlation in the R matrix. Then see if we can estimate alpha first using rescale.
    phy_ou_fn <- phytools::rescale(x = phy, model = "OU")
    likfn_ou <- function(a){
      phy_ou <- phy_ou_fn(alpha = a)
      C <- ape::vcv.phylo(phy_ou)
      ## Check if the matrix is singular
      if( matrixcalc::is.singular.matrix(x = C, tol = min_tol) ) return( Inf )
      tt <- phytools::phyl.vcv(X = x, C = C, lambda = 1)
      return( -1 * ratematrix::likelihoodFunction(data = x, phy = phy_ou, R = tt$R, root = tt$alpha) )
    }
    ## Find maximum value for alpha:
    fn_to_search_ou <- function(x){
      ## TRUE if it is bad.
      phy_ou <- phy_ou_fn(alpha = x, sigsq = 1)
      return( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_ou), tol = min_tol) )
    }
    bound <- binary_search(fn = fn_to_search_ou, lb = 0, ub = 1000, chunk = 0.000001)
    ## Estimate MLE for alpha parameter:
    fit_l <- nloptr::nloptr(x0 = bound[1]/2, eval_f = likfn_ou
                            , lb = .Machine$double.xmin, ub = bound[1]
                            , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                          , xtol_rel = 1e-04, print_level = 0))
    ## This is a vector with the model specifics
    fit_model_solution <- list(model = "OU", alpha = fit_l$solution, npar = 1
                               , loglik = -1*fit_l$objective)

    ## Use the MLE of the lambda to compute the probabilistic PCA model.
    phy_ou <- phy_ou_fn(alpha = fit_model_solution[[2]])
    ## Compute the pics with the rescaled tree:
    pic <- apply(X = x, MARGIN = 2, FUN = pic, phy = phy_ou)
  }
  if( model == "EB" ){
    phy_eb_fn <- phytools::rescale(x = phy, model = "EB")
    likfn_eb <- function(a){
      phy_eb <- phy_eb_fn(a = a)
      C <- ape::vcv.phylo(phy_eb)
      ## Need to have another check for matrix singularity:
      if( matrixcalc::is.singular.matrix(x = C, tol = min_tol) ) return( Inf )
      tt <- phytools::phyl.vcv(X = x, C = C, lambda = 1)
      return( -1 * ratematrix::likelihoodFunction(data = x, phy = phy_eb, R = tt$R
                                                  , root = tt$alpha) )
    }

    ## Search for the minimum value for EB parameter.
    fn_to_search_eb <- function(x){
      phy_eb <- phy_eb_fn(a = -1 * x) ## Searching negative space.
      ## Need to have another check for matrix singularity:
      return( matrixcalc::is.singular.matrix(x = ape::vcv.phylo(phy_eb), tol = min_tol ) )
    }
    bound <- -1*binary_search(fn = fn_to_search_eb, lb = 0, ub = 1000, chunk = 0.000001)
    ## Now do the likelihood search:
    fit_l <- nloptr::nloptr(x0 = bound[1]/2, eval_f = likfn_eb, lb = bound[1], ub = 0
                            , opts = list(algorithm = "NLOPT_LN_SBPLX"
                                          , xtol_rel = 1e-04, print_level = 0))
    ## This is a vector with the model specifics
    fit_model_solution <- list(model = "EB", a = fit_l$solution, npar = 1
                               , loglik = -1*fit_l$objective)

    ## Use the MLE of the lambda to compute the probabilistic PCA model.
    phy_eb <- phy_eb_fn(a = fit_model_solution[[2]])
    ## Compute the pics with the rescaled tree:
    pic <- apply(X = x, MARGIN = 2, FUN = pic, phy = phy_eb)
  }

  ## Now conduct probabilistic PCA with the independent contrasts:
  return( ProbPCA(x = pic, ret_dim = ret_dim, scale = scale) )
}

