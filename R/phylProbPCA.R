#' @title Probabilistic Phylogenetic PCA
#'
#' @description Function to perform probabilistic phylogenetic PCA. Allows for fit of alternative models of trait evolution using branch length transformation.
#' @param phy phylogeny in "phylo" format.
#' @param x a matrix with traits in columns and species values in rows. Rownames must match the tip labels of phylogeny.
#' @param ret_dim number of dimensions (PC axes) to be kept by the model.
#' @param model choice of model of trait evolution. One of "BM", "lambda", "OU", or "EB".
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
#' @examples
#' phy <- ratematrix::anoles$phy[[1]]
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#' ppca <- phylProbPCA(phy = phy, x = dt, ret_dim = 2)
#' doBiplot(x = ppca, add_margin = 0.3)
#'
#' @export
phylProbPCA <- function(phy, x, ret_dim = 2, model = "BM", quiet = FALSE){
  ## Check if x is a matrix object.
  if( !inherits(x = x, what = "matrix") ){
    stop("x needs to be a matrix. Cannot be a data.frame.")
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
    C <- ape::vcv.phylo(phy)
    invC <- solve(C)
    a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
    A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
    R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
    muhat <- c(a)
    fit_model_solution <- list(model = "BM")
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
    C <- ape::vcv.phylo(phy_lambda)
    invC <- solve(C)
    a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
    A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
    R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
    muhat <- c(a)
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
    C <- ape::vcv.phylo(phy_ou)
    invC <- solve(C)
    a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
    A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
    R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
    muhat <- c(a)
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
    C <- ape::vcv.phylo(phy_eb)
    invC <- solve(C)
    a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
    A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
    R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
    muhat <- c(a)
  }

  ## Adding BM parameters to the model output:
  fit_model_solution$R <- R
  fit_model_solution$muhat <- muhat

  eigR <- eigen(R, TRUE)

  mlsig2 <- (sum(eigR$values[(q+1):d])) / (d-q)
  ## Tipping and Bishop MLE formula include the square root of the eigenvalues.
  mlW <- (eigR$vectors[,1:q]) %*% (diag(sqrt(eigR$values[1:q] - mlsig2)))

  M <- (t(mlW)%*%mlW)+(diag(ncol(mlW))*mlsig2)
  ## aux.bicgstab_fork is a fork of Rdimtools:::aux.bicgstab version 1.1.2
  SOL <- aux.bicgstab_fork(M, t(mlW), verbose = FALSE)
  ## aux.adjprojection_fork is a fork of Rdimtools:::aux.adjprojection version 1.1.2
  projection <- aux.adjprojection_fork(t(SOL$x))

  ## This should be changed to the likelihood of the pppca as we describe in the paper.
  lik <- ratematrix::likelihoodFunction(data = x, phy = phy
                                        , R = mlW%*%t(mlW)+mlsig2*diag(d)
                                        , root = muhat)

  ## Number of parameters depend on model:
  if( model == "BM" ){
    params <- (d*q) - (0.5*q*(q-1)) + 1
  } else{
    ## One extra parameter for all other models.
    params <- (d*q) - (0.5*q*(q-1)) + 2
  }

  ## Save the MLE parameters:
  par_list <- list(R = mlW%*%t(mlW)+mlsig2*diag(d), root = muhat, npar = params)

  AIC <- (-2*lik) + (2*params)
  AICc <- (-2*lik) + (2*params*(N/(N-params-1)))
  BIC <- (-2*lik) + (params*log(N))
  # DANIEL : 7. need to compute the scores given the projection
  ## The score computation will need to change for the phylogenetic PCA.
  Xc <- sweep(x, 2, muhat, "-")
  scores <- t( t(projection) %*% t(Xc) )
  colnames(scores) <- paste0("PC", 1:q)

  ## Compute the loadings:
  Ccv <- t(x-A) %*% invC %*% scores / (N-1)
  L <- matrix(data = NA, nrow = d, ncol = q)
  colnames(L) <- paste("PC",1:q,sep="")
  rownames(L) <- colnames(x)
  for(i in 1:d){
    for(j in 1:q){
      L[i,j] <- Ccv[i,j] / sqrt( R[i,i] * eigR$values[j] )
    }
  }

  ## Computation of eigenvalues given projection (eigenvectors).
  e_val <- vector(mode = "numeric", length = ncol(projection))
  for( i in 1:length(e_val) ){
    lm_dt <- data.frame(R_proj = R %*% projection[,i], proj = projection[,i])
    lm_tmp <- stats::lm(R_proj~0+proj, data = lm_dt)
    e_val[i] <- unname(stats::coef(lm_tmp))
  }

  ## RETURN
  ## g is always 1, unless it is a mixed model
  result <- list(scores = scores, e_values = e_val, e_vectors = projection
                 , model_fit = fit_model_solution, loadings = L
                 , varnames = colnames(x), sig = mlsig2, mle.W = mlW, AIC = AIC
                 , AICc = AICc, BIC = BIC, loglik = lik, g = 1, mle_par_list = par_list)
  class(result) <- c(class(result), "phylPPCA")
  return(result)
}
