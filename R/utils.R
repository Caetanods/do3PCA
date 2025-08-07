## Test
## Dependency of some unexported functions from Rdimtools:

#' @noRd
#' @import Rdimtools
#' @importFrom utils getFromNamespace
linsolve_Rdimtools <- utils::getFromNamespace(x = "linsolve.bicgstab.single", ns = "Rdimtools")

#' @noRd
#' @import Rdimtools
#' @importFrom utils getFromNamespace
linsolve_sparse_Rdimtools <- utils::getFromNamespace(x = "linsolve.bicgstab.single.sparse", ns = "Rdimtools")

#' @noRd
#' @import Rdimtools
#' @importFrom utils getFromNamespace
aux.adjprojection_fork <- utils::getFromNamespace(x = "aux.adjprojection", ns = "Rdimtools")

#' @noRd
#' @import Rdimtools
#' @importFrom utils getFromNamespace
aux.bicgstab_fork <- utils::getFromNamespace(x = "aux.bicgstab", ns = "Rdimtools")

binary_search <- function(fn, lb, ub, chunk = 0.00001, max_reps = 100000){
  ## Use binary search to find an interval of length chunk in which lb is FALSE and ub is TRUE.
  ## Assumes fn will return TRUE or FALSE and that TRUE values should be within high numbers.
  i <- 0
  while(i <= max_reps){
    i <- i+1
    center <- mean(c(ub, lb))
    ## We want inter[1] finite and inter[2] Inf.
    inter <- c(center-(chunk/2), center+(chunk/2))
    if( !fn(inter[1]) ){
      ## Check if finished:
      if( fn(inter[2]) ) break()
      ## If not, choose right side:
      lb <- center
    } else{
      ## Choose left side:
      ub <- center
    }
  }
  if( i > max_reps ) warning("Max reps reached, result might be unstable.")
  return( c(lb,ub) )
}

#' @noRd
#' @importFrom ape vcv.phylo
#' @importFrom Matrix solve
fastR <- function(x, phy){
  C <- ape::vcv.phylo(phy)
  invC <- Matrix::solve(C, sparse = TRUE)
  a <- matrix(colSums(invC %*% x) / sum(invC), ncol(x), 1)
  A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
  R <- t(x-A) %*% invC %*% (x-A) / (nrow(C)-1)
}

#' @noRd
#' @importFrom ape vcv.phylo
#' @importFrom stats lm coef
getPhyloProjection <- function(x, invC, R, R_decom, W, a, N, sigma, d, q, squared = FALSE){
  # x = data
  # phy = phylogenetic tree.
  # W = the W matrix.
  # a = the phylogenetic mean.
  # sigma = can be squared or not
  # squared = indicate if sigma is squared or not

  if( !squared ){
    sigma <- sigma^2
  }
  ## Following the projection computed in do3PCA::phylProbPCA.
  M <- (t(W)%*%W)+(diag(ncol(W))*sigma)
  SOL <- aux.bicgstab_fork(M, t(W), verbose=FALSE)
  projection <- aux.adjprojection_fork(t(SOL$x))
  Xc <- sweep(x, 2, a, "-")
  scores <- t( t(projection) %*% t(Xc) )

  ## Compute the loadings:
  A <- matrix(rep(a, nrow(x)), nrow(x), ncol(x), byrow=T)
  # invC <- solve(ape::vcv.phylo(phy))
  Ccv <- t(x-A) %*% invC %*% scores / (N-1)
  Lm <- matrix(data = NA, nrow = d, ncol = q)
  colnames(Lm) <- paste("PC",1:q,sep="")
  rownames(Lm) <- colnames(x)
  for(i in 1:d){
    for(j in 1:q){
      Lm[i,j] <- Ccv[i,j] / sqrt( R[i,i] * R_decom$values[j] )
    }
  }

  ## Computation of eigenvalues given projection (eigenvectors).
  e_val <- vector(mode = "numeric", length = ncol(projection))
  for( i in 1:length(e_val) ){
    lm_dt <- data.frame(R_proj = R %*% projection[,i], proj = projection[,i])
    lm_tmp <- stats::lm(R_proj~0+proj, data = lm_dt)
    e_val[i] <- unname(stats::coef(lm_tmp))
  }

  result <- list()
  result$scores <- scores
  result$e_values <- e_val
  result$projection <- projection
  result$sig <- sigma
  result$mle.W <- W
  result$varnames <- colnames(x)
  return( result )
}

#' @noRd
#' @importFrom stats lm coef
getProjection <- function(x, W, S, sigma, squared = FALSE){
  # x = data
  # W = the W matrix.
  # sigma = can be squared or not
  # squared = indicate if sigma is squared or not
  if( !squared ){
    sigma <- sigma^2
  }
  M <- (t(W)%*%W)+(diag(ncol(W))*sigma)
  SOL <- aux.bicgstab_fork(M, t(W), verbose=FALSE)
  projection <- aux.adjprojection_fork(t(SOL$x))
  ## Find the eigenvalues given W:
  e_val <- vector(mode = "numeric", length = ncol(projection))
  for( i in 1:length(e_val) ){
    lm_dt <- data.frame(S_proj = S %*% projection[,i], proj = projection[,i])
    lm_tmp <- stats::lm(S_proj~0+proj, data = lm_dt)
    e_val[i] <- unname(stats::coef(lm_tmp))
  }
  result <- list()
  result$scores <- (x%*%projection)
  result$projection <- projection
  result$e_values <- e_val
  result$sig <- sigma
  result$mle.W <- W
  result$varnames <- colnames(x)
  return( result )
}
