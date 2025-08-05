#' Broken stick test for number of principal components
#'
#' Function generates the expected distribution of variance across the principal components and compares it with the observed vector of explained variance (if provided).
#' @param n_dim the number of variables of the data (required).
#' @param var_vec the vector of explained variance of the PCA (optional).
#' @param tol a tolerance value that is used to check if the observed explained variance is larger than the expected under the Broken Stick model.
#' @param plot if a plot should be produced.
#'
#' @returns Function will return the vector of explained variance for each PC. Will also return the number of PCs to keep and optionally a plot of the results if 'var_vec' is provided.
#' @importFrom graphics barplot lines
#' @export
BrokenStick <- function(n_dim, var_vec = NULL, tol = 0.0, plot = TRUE){
  ## Add this helping function to the package.
  x <- 1/(1:n_dim)
  bs_vec <- sapply(1:n_dim, function(k) sum(x[k:n_dim])/n_dim)
  if( !is.null(var_vec) ){
    if( !length(var_vec) == n_dim ){
      stop("length of var_vec must equal n_dim.")
    }
    std_var <- var_vec / sum(var_vec)
    if( plot ){
      ## Make a plot of the results:
      midbar <- barplot(height = bs_vec, names.arg = 1:n_dim, col = "lightgrey"
                        , border = FALSE, ylim = c(0,1))
      lines(x = midbar, y = std_var, type = "b", col = "red", lwd = 2)
    }
    test_pcs <- ((std_var - bs_vec) - tol) > 0.0
    test_TRUE <- TRUE
    i <- 0
    while( test_TRUE ){
      i <- i+1
      test_TRUE <- test_pcs[i]
    }
    ## First FALSE position is the number of returned PCs:
    return( list(broken_stick = bs_vec, n_pc_to_keep = i-1) )
  } else{
    return( list(broken_stick = bs_vec) )
  }
}
