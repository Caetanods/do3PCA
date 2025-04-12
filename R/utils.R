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
