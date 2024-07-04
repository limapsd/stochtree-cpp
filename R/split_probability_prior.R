
#' Title
#'
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
log_sum_exp = function(v){
  
  n = length(v)
  mx = max(v)
  sm = 0
  for(i in 1:n){
    sm = sm + exp(v[i] - mx)
  }
  return(mx + log(sm))
}




log_gamma_sampler = function(shape){
  
  y = log(rgamma(1, shape+ 1))
  z = log(runif(1))/shape
  return(y+z)
}

#' Title
#'
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
log_dirichilet_sampler = function(alpha){
  
  k = length(alpha)
  draw = rep(0,k)
  for(j in 1:k){
    draw[j] = log_gamma_sampler(alpha[j])
  }
  lse = log_sum_exp(draw)
  for(j in 1:k){
    draw[j] = draw[j] - lse
  }
  return(draw)
}


#' Title
#'
#' @param nv Vector with the number of splits for each variable
#' @param m_i Endogenous variable index
#' @param lag_index index for choosing the appropiate weight in the dirichilet
#' @param sigma_ols calibrated weight for the cross-lags
#' @param lambda_1 shrinkage parameter for i = j
#' @param lambda_2 shrinkage parameter for i != j
#'
#' @return List of probabilities according to the Minessota criteria.
#' @export
#'
#' @examples
draw_minessota_split = function(nv, m_i, lag_index, sigma_ols, lambda_1, lambda_2){
  
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
    lpv = log_dirichilet_sampler(phi_)
    return(lpv)

}




#' Title
#'
#' @param nv 
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
draw_dart_splits = function(nv,theta = 1){
  n   = length(nv)
  theta_ = rep(0, n)
  for(i in 1:n){
    theta_[i] = theta/n + nv[i]
  }
  lpv = log_dirichilet_sampler(theta_)
  return(lpv)
}

discrete = function(wts) {
  p <- length(wts)
  x <- 0
  vOut <- rep(0, p)
  vOut <- rmultinom(1, size = 1, prob = wts)
  if (vOut[1] == 0) {
    for (j in 2:p) {
      x <- x + j * vOut[j]
    }
  }
  return(x)
}

#' Title
#'
#' @param theta 
#' @param lpv 
#' @param a 
#' @param b 
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
draw_theta_update = function(theta, lpv, a , b, rho) {
  
  p      = length(lpv)
  sumlpv = sum(lpv)
  lambda_g <- seq(1 / 1001, 1000 / 1001, length.out = 1000)
  theta_g <- lambda_g * rho / (1 - lambda_g)
  lwt_g    = rep(0, 1000)
  
  for (k in 1:1000) {
    theta_log_lik = lgamma(theta_g[k]) - p * lgamma(theta_g[k] / p) + (theta_g[k] / p) * sumlpv
    beta_log_prior = (a - 1) * log(lambda_g[k]) + (b - 1) * log(1 - lambda_g[k])
    lwt_g[k] = theta_log_lik + beta_log_prior
  }
  
  lse <- log_sum_exp(lwt_g)
  lwt_g <- exp(lwt_g - lse)
  weights <- lwt_g / sum(lwt_g)
  theta <- theta_g[discrete(weights)]
  
  return(theta)
}

