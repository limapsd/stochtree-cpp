

# log_sum_exp <- function(v){
#   
#   n = length(v)
#   mx = max(v)
#   sm = 0
#   for(i in 1:n){
#     sm = sm + exp(v[i] - mx)
#   }
#   return(mx + log(sm))
# }

#' Title
#'
#' @param v 
#'
#' @return
#' @export
#'
#' @examples
log_sum_exp = function(v){
  
  M = max(v)
  r = M + log(sum(exp(v - M)))
  return(r)
}

rho_to_alpha <- function(rho, scale){
  return ( (scale* rho)/(1.0 -rho))
}

alpha_to_rho <- function(alpha, scale){
  return (alpha/(alpha + scale))
}

logpdf_beta <- function(x, a, b){
  r = (a-1.0) * log(x) + (b-1.0) * log(1 - x) - lbeta(a,b)
  return (r)
}

rho_loglik <- function(mean_log_s, p, alpha_scale, alpha_shape_1, alpha_shape_2){
  
  loglike = function(rho){
    alpha = rho_to_alpha(rho, alpha_scale)
    r = alpha * mean_log_s + lgamma(alpha) - p*lgamma(alpha/p) + logpdf_beta(rho, alpha_shape_1, alpha_shape_2)
    return(r)
  }
  return(loglike)
}


sample_class <- function(probs) {
  U <- runif(1)  # Generate a single uniform random number between 0 and 1
  foo <- 0.0
  K <- length(probs)  # Length of the probs vector
  
  # Sample
  for (k in 1:K) {
    foo <- foo + probs[k]
    if (U < foo) {
      return(k)
    }
  }
  return (K) # Return K if U >= sum(probs)
}



#' Title
#'
#' @param alpha 
#' @param lpv 
#' @param a 
#' @param b 
#' @param rho 
#'
#' @return
#' @export
#'
#' @examples
draw_alpha_update = function(alpha, lpv, a, b, rho){
  
  grid_size    = 1000
  rho_proposed = rep(0, grid_size - 1)
  for(i in 1:(grid_size - 1)){
    rho_proposed[i] = i/grid_size
  }
  loglikes = rep(0,grid_size-1)

  p = length(lpv)
  
  alpha_scale   =  rho
  alpha_shape_1 =  a
  alpha_shape_2 =  b

  loglike  = rho_loglik(mean(lpv), p, alpha_scale, alpha_shape_1, alpha_shape_2)
  #grid_size-1
  for(i in 1:grid_size-1){
    loglikes[i] = loglike(rho_proposed[i])
  }
 
  loglikes = exp(loglikes - log_sum_exp(loglikes)) 
  
  # index = sample_class(loglikes)
  index <- sample(1:length(loglikes), size = 1, prob = loglikes, replace = F)
  rho_up = rho_proposed[index]
  alpha = rho_to_alpha(rho_up, alpha_scale)
  return(alpha)

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
draw_theta_update <- function(theta, lpv, a , b, rho) {

  p      = length(lpv)
  sumlpv = sum(lpv)

  grid_size = 100
  lambda_g  = rep(0, grid_size)
  theta_g   = rep(0, grid_size)
  lwt_g     = rep(0, grid_size)

  for(k in 1:(grid_size)){
    lambda_g[k] = k/(grid_size +1)
    theta_g[k] = (lambda_g[k]*rho)/(1.-lambda_g[k])
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



log_gamma_sampler <- function(shape) {

  if (shape >= 0.1){
    # gamma_sample = rgamma(1, shape, 1)
    gamma_sample = rgamma(1, shape +1.0, 1)
    r = log(gamma_sample) + log(runif(1.0))/shape
    # r = log(gamma_sample)
    return(r)
    }

  a = shape
  L = 1.0/a - 1.0
  w = exp(-1.0) * a / (1.0 - a)
  ww = 1.0 / (1.0 + w)
  z = 0.0
  repeat {
    U = runif(1)
    if (U <= ww) {
      z = -log(U / ww)
    } else {
      z = log(runif(1)) / L
    }
    eta = ifelse(z >= 0, -z, log(w) + log(L) + L * z)
    h = -z - exp(-z / a)
    if (h - eta > log(runif(1))) break
  }
  # cat("Shape: ", shape, "\n")
  # cat("Sample: ", -z/a, "\n")

  return(-z/a)
}


# 
# log_gamma_sampler = function(shape){
#   y = log(rgamma(1, shape+ 1, 1))
#   z = log(runif(1))/shape
#   return(y+z)
# }

#' Title
#'
#' @param alpha 
#'
#' @return
#' @export
#'
#' @examples
log_dirichilet_sampler <- function(alpha){

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




#' Title
#'
#' @param nv 
#' @param theta 
#'
#' @return
#' @export
#'
#' @examples
draw_dart_splits <- function(nv,theta = 1){
  n   = length(nv)
  theta_ = rep(0, n)
  for(i in 1:n){
    theta_[i] = theta/n + nv[i]
  }
  # theta_ = theta/n + nv
  lpv = log_dirichilet_sampler(theta_)
  return(lpv)
}

discrete <- function(wts) {
  p <- length(wts)
  x <- 0
  vOut <- rep(0, p)
  vOut <- rmultinom(1, size = 1, prob = wts)
  if (vOut[1] == 0) {
    for (j in 2:p) {
      x <- x + (j-1) * vOut[j]
    }
  }
  return(x)
}

