#' Sample one iteration of the split probabilities following Linero(2018) 
#'
#' @param variable_count_splits number of splits for each variable 
#' @param alpha scale of the dirichilet distribtuion
#' @param rng C++ random number generator
#'
#' @export
sample_dart_splits_one_iteration <- function(variable_count_splits,alpha, rng) {
  results = sample_dart_splits_one_iteration_cpp(variable_count_splits, alpha, rng$rng_ptr)
  r = list("lse_dir" =  results[[1]], "lpv" =  results[[2]]) 
  
  return(r)
}


#'  Sample one iteration of the scale of the Dirichlet Distribution
#'
#' @param log_prob_vector 
#' @param alpha 
#' @param a 
#' @param b 
#' @param rho 
#' @param rng C++ random number generator
#'
#' @return
#' @export
#'
#' @examples
sample_alpha_one_iteration <- function(log_prob_vector, a , b, rho, rng){
    result = sample_alpha_one_iteration_cpp(log_prob_vector, a, b, rho, rng$rng_ptr)
    r = list("loglikes" = result[[1]], "alpha" = result[[2]], "lse_alpha" = result[[3]])
    
    return(r)
}



#'  Sample one iteration of the Minnesota Dart Proposal
#'
#' @param phi 
#' @param rng 
#'
#' @return
#' @export
#'
#' @examples
sample_minnesota_one_iteration <- function(phi,rng){
  result = sample_minnesota_dart_one_iteration_cpp(phi, rng$rng_ptr)
  r = list("lse_dir" =  results[[1]], "lpv" =  results[[2]]) 
  return(r)
}

#' Title
#'
#' @param nv Vector with the number of splits for each variable
#' @param m_i Endogenous variable index
#' @param lag_index index for choosing the appropiate weight in the dirichilet
#' @param sigma_ols calibrated weight for the cross-lags
#' @param lambda_1 shrinkage parameter for i = j
#' @param lambda_2 shrinkage parameter for i != j
#' @param rng 
#' @return List of probabilities according to the Minessota criteria.
#' @export
#'
#' @examples
draw_minessota_split <- function(nv, m_i, lag_index, sigma_ols, lambda_1, lambda_2, rng){
  
    n = length(nv)
    phi_ = rep(0, n)
    for(i in 1:n){
        if(any(i == lag_index[m_i,])){
            l = which(lag_index[m_i,] == i)
            phi_[i] = lambda_1/(l^2) + nv[i]
        }else{
            var_search = which(lag_index == i ,arr.ind=TRUE)
            l  = var_search[2]
            m_j =  lag_index[var_search[1],1]
            phi_[i] = (lambda_2 * sigma_ols[m_i])/(l^2 * sigma_ols[m_j]) + nv[i]
            
        }
    }
    # lpv = log_dirichilet_sampler(phi_)
     minn_dart = sample_minnesota_one_iteration(phi_,rng)
    return(minn_dart$lpv)

}

