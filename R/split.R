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
