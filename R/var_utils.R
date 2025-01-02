
# lag variables

#' Title
#'
#' @param X 
#' @param lag 
#'
#' @return
#' @export
#'
#' @examples
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}



                      
  


preprocessmBartParams <- function(params) {

    # Default parameter values
  processed_params <- list(alpha = 0.95, beta = 2.0, min_samples_leaf = 5,
                      max_depth = 10, n_trees = 200, cutpoint_grid_size = 100,
                      nu = 4, lambda = 0.5, a_leaf = 2, b_leaf = NULL, 
                      a_forest = 1, b_forest = 1,random_seed = -1,sample_sigma_leaf = TRUE, 
                      alpha_dart = NULL, rho_dart = NULL, a_dart = 1, b_dart = 0.5)
  
  for (key in names(params)) {
    if (!key %in% names(processed_params)) {
      stop ("Variable", key, "is not valid BART parameter")
    }
    val <-params[[key]]
    if (!is.null(val)) processed_params[[key]] <- val
  }

  return(processed_params)
  
}
                      