## Function to get a posterior distribution of scores to plot.
getPostProjection <- function(post, R_decom, q, x, invC, R, W, a, N, d, post_n){
  ## post is in log:
  post <- exp(post)
  ## Take the samples from the posterior distribution:
  sample_id <- sample(x = 1:nrow(post), size = post_n)
  post <- post[sample_id,]
  
  ## Make the function to compute the posterior scores:
  fn_factory <- function(R_decom, q, x, invC, R, W, a, N, d) function(L, sigma){
    ## Reconstruction of W matrix from L vector:
    ## With the W matrix and sigma we can get the principal components.
    L_q <- diag(x = L)
    U_q <- R_decom$vectors[,1:q]
    I_q <- diag(nrow = q, ncol = q)
    W <- U_q %*% (L_q - (sigma^2*I_q))^0.5 %*% I_q
    
    ## Compute the projection:
    mcmc_proj <- getPhyloProjection(x = x, invC = invC, R = R, R_decom = R_decom
                                    , W = W, a = a, N = N, d = d, q = q
                                    , sigma = sigma)
    
    ## Compute the percent explained variance of the axes:
    mcmc_pc_var <- mcmc_proj$e_values / sum(mcmc_proj$e_values) * 100
    ## Compute the loadings:
    mcmc_loadings <- mcmc_proj$projection^2 * 100
    rownames(mcmc_loadings) <- mcmc_proj$varnames
    colnames(mcmc_loadings) <- paste0("loadings_PC", 1:q)
    ## Get the scores:
    scores <- mcmc_proj$scores
    colnames(scores) <- paste0("PC_", 1:q)
    
    return( list(scores = scores, pc_var = mcmc_pc_var
                 , loadings = mcmc_loadings) )
  }
  
  ## Generate the function with the constants:
  proj_fn <- fn_factory(R_decom = R_decom, q = q, x = x, invC = invC, R = R
                        , W = W, a = a, N = N, d = d)
  ## Loop over the posterior samples:
  out_list <- list()
  for( i in 1:nrow(post) ){
    out_list[[i]] <- proj_fn(L = post[i,1:(ncol(post)-1)], sigma = post[i,ncol(post)])
  }
  
  ## Return the list of scores:
  return( out_list )
}
