#' @title Make biplot for any type of PCA
#'
#' @description Function to make biplots for any kind of PCA. It accepts the outputs from standard PCA (princomp and prcomp). It also works with the "phylProbPCA" and "ProbPCA" functions. It provides more options to the plot than the standard "stats::biplot".
#' @param x output from PCA analysis.
#' @param choices numeric vector of length 2. Use to choose which of the PC axes to plot. Default plots first and second axes: "choices = c(1,2)".
#' @param scale numeric value between 0 and 1. Same as in "stats::biplot.princomp". See ?biplot.princomp for more information.
#' @param pc.biplot logical. If TRUE it will produce a "principal component biplot" (sensu Gabriel, 1971). Same as in "stats::biplot.princomp". See ?biplot.princomp for more information.
#' @param col character vector of length 3 with the colors of the biplot. First color is used for the score points (or sample sames), second color for arrows and variable names, and third color for the right and top-side ticks (plot axes).
#' @param ... extra parameters for the function. Same as "stats::biplot".
#' @returns makes a biplot of the PCA results.
#' @details
#' Function has the same options as "stats::biplot", with the addition of the following arguments. "plot_dimnames" controls is the names of the samples (species) will be plotted. "add_points" controls if the score points will be plotted. "add_margin" is a numeric value that expands the area of the plot. You can use this to make sure the names of variables and samples (species) fit the plot.
#'
#' @examples
#' phy <- ratematrix::anoles$phy[[1]]
#' dt <- as.matrix( ratematrix::anoles$data[,1:3] )
#'
#' ## Using probabilistic phylogenetic PCA
#' phylppca <- phylProbPCA(phy = phy, x = dt, ret_dim = 2)
#' doBiplot(x = phylppca, add_margin = 0.3)
#'
#' ## Using standard phylogenetic PCA
#' phylpca <- phytools::phyl.pca(tree = phy, Y = dt)
#' doBiplot(x = phylpca, add_margin = 0.3)
#'
#' ## Using probabilistic PCA
#' ppca <- ProbPCA(x = dt)
#' doBiplot(x = ppca, add_margin = 0.3)
#'
#' ## Using standard PCA
#' pca1 <- princomp(x = dt)
#' doBiplot(x = pca1, add_margin = 0.1)
#'
#' ## Using standard PCA
#' pca2 <- prcomp(x = dt)
#' doBiplot(x = pca2, add_margin = 0.1)
#'
#' @export
doBiplot <- function(x, choices = 1L:2L, scale = 1, pc.biplot = FALSE, col, ...){
  ## Wrapper function to do biplots for each of the types of PCA.
  ## This could be a S3 method, but it is not.
  if( inherits(x = x, what = "princomp") ){
    biplot_princomp(x, choices = choices, scale = scale, pc.biplot = pc.biplot, col = col, ...)
  } else if( inherits(x = x, what = "prcomp") ){
    biplot_prcomp(x, choices = choices, scale = scale, pc.biplot = pc.biplot, col = col, ...)
  } else if( inherits(x = x, what = "phylPPCA") ){
    biplot_phylo_ppca(x, choices = choices, scale = scale, pc.biplot = pc.biplot, col = col, ...)
  } else if( inherits(x = x, what = "PPCA") ){
    biplot_ppca(x, choices = choices, scale = scale, pc.biplot = pc.biplot, col = col, ...)
  } else if( inherits(x = x, what = "phyl.pca") ){
    biplot_phylo_pca(x, choices = choices, scale = scale, pc.biplot = pc.biplot, col = col, ...)
  } else{
    stop("Function works with object class of type 'princomp', 'prcomp', 'phylPPCA', 'PPCA', or 'phyl.pca'.")
  }
}

#' @importFrom graphics par text points axis box arrows
#' @importFrom grDevices dev.hold dev.flush
#' @importFrom stats sd
#' @noRd
biplot_default <-
  function(x, y, var.axes = TRUE, col, cex = rep(graphics::par("cex"), 2),
           xlabs = NULL, ylabs = NULL, expand=1, xlim = NULL, ylim = NULL,
           arrow.len = 0.1, main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
           plot_dimnames = FALSE, add_points = TRUE, add_margin = 0, ...)
  {
    n <- nrow(x)
    p <- nrow(y)
    if(missing(xlabs)) {
      xlabs <- dimnames(x)[[1L]]
      if(is.null(xlabs)) xlabs <- 1L:n
    }
    xlabs <- as.character(xlabs)
    dimnames(x) <- list(xlabs, dimnames(x)[[2L]])
    if(missing(ylabs)) {
      ylabs <- dimnames(y)[[1L]]
      if(is.null(ylabs)) ylabs <- paste("Var", 1L:p)
    }
    ylabs <- as.character(ylabs)
    dimnames(y) <- list(ylabs, dimnames(y)[[2L]])

    if(length(cex) == 1L) cex <- c(cex, cex)
    # if(missing(col)) {
    #   col <- par("col")
    #   if (!is.numeric(col)) col <- match(col, palette(), nomatch=1L)
    #   col <- c(col, col + 1L, col + 1L)
    # }
    # else if(length(col) == 1L) col <- c(col, col)
    if(missing(col)) {
      ## Expanding the col types.
      col <- graphics::par("col")
      col <- c("grey", "#DF536B", col)
    }
    else if(length(col) != 3L){
      warning( "Using default colors. col argument needs to be a vector with 3 elements." )
      col <- graphics::par("col")
      col <- c("grey", "#DF536B", col)
    }

    unsigned.range <- function(x)
      c(-abs(min(x, na.rm=TRUE)), abs(max(x, na.rm=TRUE)))
    rangx1 <- unsigned.range(x[, 1L])
    rangx2 <- unsigned.range(x[, 2L])
    rangy1 <- unsigned.range(y[, 1L])
    rangy2 <- unsigned.range(y[, 2L])

    if(missing(xlim) && missing(ylim))
      xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
    else if(missing(xlim)) xlim <- rangx1
    else if(missing(ylim)) ylim <- rangx2
    ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
    ## Setting graphics parameters to reset at function exit:
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(oldpar))
    ## on.exit(graphics::par(op))
    ## op <- graphics::par(pty = "s")
    ## Add to margins before starting to plot:
    xlim <- c(xlim[1] - add_margin, xlim[2] + add_margin)
    ylim <- c(ylim[1] - add_margin, ylim[2] + add_margin)
    ## if(!is.null(main))
      ## op <- c(op, graphics::par(mar = graphics::par("mar")+c(0,0,1,0)))
    plot(x, type = "n", xlim = xlim, ylim = ylim, col = col[1L],
         xlab = xlab, ylab = ylab, sub = sub, main = main, ...)
    if( plot_dimnames ){ ## Allow taking out dimnames
      graphics::text(x, xlabs, cex = cex[1L], col = col[1L], ...)
    }
    if( add_points ){ ## Can add the points in the plot.
      graphics::points(x, pch = 19, col = col[1L])
    }
    graphics::par(new = TRUE)
    grDevices::dev.hold(); on.exit(grDevices::dev.flush(), add = TRUE)
    plot(y, axes = FALSE, type = "n", xlim = xlim*ratio, ylim = ylim*ratio,
         xlab = "", ylab = "", col = col[1L], ...)
    # axis(3, col = col[2L], ...)
    # axis(4, col = col[2L], ...)
    graphics::axis(3, col = col[3L], ...)
    graphics::axis(4, col = col[3L], ...)
    graphics::box(col = col[1L])
    graphics::text(y, labels=ylabs, cex = cex[2L], col = col[2L], ...)
    if(var.axes)
      graphics::arrows(0, 0, y[,1L] * 0.8, y[,2L] * 0.8, col = col[2L], length=arrow.len)
    invisible()
}

## For the princomp function
biplot_princomp <- function(x, choices = 1L:2L, scale = 1, pc.biplot=FALSE, col, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  lam <- x$sdev[choices]
  if(is.null(n <- x$n.obs)) n <- 1
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  biplot_default(t(t(scores[, choices]) / lam),
                 t(t(x$loadings[, choices]) * lam), col = col, ...)
  invisible()
}

## For the prcomp function
biplot_prcomp <- function(x, choices = 1L:2L, scale = 1, pc.biplot=FALSE, col, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$x))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  if(is.complex(scores))
    stop("biplots are not defined for complex PCA")
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  biplot_default(t(t(scores[, choices]) / lam),
                 t(t(x$rotation[, choices]) * lam), col = col, ...)
  invisible()
}

## Works with Rdimtools format.

#' @importFrom stats sd
#' @noRd
biplot_ppca <- function(x, choices = 1L:2L, scale = 1, pc.biplot = FALSE, col, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  ## Problem in the next line. The scores can be a list! Only matrix if there is a single group.
  if(is.complex(scores))
    stop("biplots are not defined for complex PCA")
  lam <- apply(scores, 2, stats::sd)
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  ## DANIEL: Adding step to prevent wrong operations.
  lam_matrix <- matrix(lam, nrow = n, ncol = length(lam), byrow = TRUE)
  ## DANIEL: For some reason the scores matrix is in the incorrect orientation.
  if( ncol(scores) == ncol(lam_matrix) ){
    ## This should be always the option.
    x_plot <- scores / lam_matrix
  } else{
    warning( "unmatched ncol of scores")
    x_plot <- t(t(scores) / lam_matrix)
  }
  n_var <- nrow(x$projection)
  lam_matrix <- matrix(lam, nrow = n_var, ncol = length(lam), byrow = TRUE)
  loadings_mat <- x$projection
  if( ncol(loadings_mat) == ncol(lam_matrix) ){
    y_plot <- loadings_mat * lam_matrix
  } else{
    y_plot <- t(t(loadings_mat) * lam_matrix)
  }
  rownames(y_plot) <- x$varnames
  biplot_default(x_plot, y_plot, col = col, ...)
  invisible()
}

## This is the function for the phylogenetic probabilistic PCA.
biplot_phylo_ppca <- function(x, choices = 1L:2L, scale = 1, pc.biplot = FALSE, col, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$scores))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  ## Problem in the next line. The scores can be a list! Only matrix if there is a single group.
  if(is.complex(scores))
    stop("biplots are not defined for complex PCA")
  ## The value of "lam" is the problem!
  ## lam <- sqrt( diag(x$Eval) )
  ## lam <- apply(scores, 2, sd)
  lam <- sqrt( x$e_values )
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  ## DANIEL: Adding step to prevent wrong operations.
  lam_matrix <- matrix(lam, nrow = n, ncol = length(lam), byrow = TRUE)
  ## DANIEL: For some reason the scores matrix is in the incorrect orientation.
  if( ncol(scores) == ncol(lam_matrix) ){
    ## This should be always the option.
    x_plot <- scores / lam_matrix
  } else{
    warning( "unmatched ncol of scores")
    x_plot <- t(t(scores) / lam_matrix)
  }
  n_var <- nrow(x$e_vectors)
  lam_matrix <- matrix(lam, nrow = n_var, ncol = length(lam), byrow = TRUE)
  loadings_mat <- x$e_vectors
  if( ncol(loadings_mat) == ncol(lam_matrix) ){
    y_plot <- loadings_mat * lam_matrix
  } else{
    y_plot <- t(t(loadings_mat) * lam_matrix)
  }
  rownames(y_plot) <- x$varnames
  biplot_default(x_plot, y_plot, col = col, ...)
  invisible()
}

## Function for the phylogenetic PCA using contrasts.
biplot_phylo_pca <- function(x, choices = 1L:2L, scale = 1, pc.biplot=FALSE, col, ...)
{
  if(length(choices) != 2L) stop("length of choices must be 2")
  if(!length(scores <- x$S))
    stop(gettextf("object '%s' has no scores", deparse(substitute(x))),
         domain = NA)
  if(is.complex(scores))
    stop("biplots are not defined for complex PCA")
  lam <- sqrt( diag(x$Eval) )
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if(scale < 0 || scale > 1) warning("'scale' is outside [0, 1]")
  if(scale != 0) lam <- lam^scale else lam <- 1
  if(pc.biplot) lam <- lam / sqrt(n)
  x_plot <- t(t(scores) / lam)
  y_plot <- t(t(x$Evec) * lam)
  rownames(y_plot) <- rownames(x$Evec)
  biplot_default(x_plot, y_plot, col = col, ...)
  invisible()
}
