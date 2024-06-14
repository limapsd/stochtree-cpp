#' Title
#'
#' @param y_train 
#' @param X_train 
#' @param sv 
#' @param W_train 
#' @param group_ids_train 
#' @param rfx_basis_train 
#' @param X_test 
#' @param W_test 
#' @param group_ids_test 
#' @param rfx_basis_test 
#' @param ordered_cat_vars 
#' @param unordered_cat_vars 
#' @param cutpoint_grid_size 
#' @param tau_init 
#' @param alpha 
#' @param beta 
#' @param min_samples_leaf 
#' @param output_dimension 
#' @param is_leaf_constant 
#' @param leaf_regression 
#' @param nu 
#' @param lambda 
#' @param a_leaf 
#' @param b_leaf 
#' @param q 
#' @param sigma2_init 
#' @param num_trees 
#' @param num_gfr 
#' @param num_burnin 
#' @param num_mcmc 
#' @param random_seed 
#' @param keep_burnin 
#'
#' @return
#' @export
#'
#' @examples
varbart =  function(y_train, X_train, sv = FALSE , W_train = NULL, group_ids_train = NULL, 
                 rfx_basis_train = NULL, X_test = NULL, W_test = NULL, 
                 group_ids_test = NULL, rfx_basis_test = NULL, 
                 ordered_cat_vars = NULL, unordered_cat_vars = NULL, 
                 cutpoint_grid_size = 100, tau_init = NULL, alpha = 0.95, 
                 beta = 2.0, min_samples_leaf = 5,  output_dimension = 1,
                 is_leaf_constant = T,leaf_regression = F, nu = 3, 
                 lambda = NULL, a_leaf = 3, b_leaf = NULL, q = 0.9,
                 sigma2_init = NULL, num_trees = 250, num_gfr = 0, 
                 num_burnin = 0, num_mcmc = 100, random_seed = -1,
                 keep_burnin = F){
  
  num_samples = num_burnin + num_mcmc 
  
  M = ncol(y_train)
  TT = nrow(y_train)
  K  = ncol(X_train)
#   colnames(X) = paste(rep(colnames(Yraw), p),
#                       sort(rep(paste(".t-", sprintf("%02d", 1:p), sep = ""), each = 2)), sep = "" )
  #  Calibrate priors for sigma^2 and tau
  
  beta_ols <- solve(crossprod(X)) %*% crossprod(X,Y)
  sigma2hat = crossprod(Y-X%*%beta_ols)/TT
  quantile_cutoff <- 0.9
  if (is.null(lambda)) {
    lambda <- (sigma2hat*qgamma(1-quantile_cutoff,nu))/nu
  }
  if (is.null(sigma2_init)) sigma2_init <- sigma2hat
  if (is.null(b_leaf)) b_leaf <- var(Y)/(2*num_trees)
  if (is.null(tau_init)) tau_init <- var(Y)/(num_trees)
  current_leaf_scale <- as.matrix(tau_init)
  current_sigma2 <- sigma2_init
  
  # Data 
  forest_dataset_train_ls = c()
  
  # Random number generator (std::mt19937)
  if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
  rng <- createRNG(random_seed)
  
  # Sampling data structures
  feature_types <- as.integer(feature_types)
  #outcome_train = createOutcome(Y[,1])
  
  # Container of forest samples
  forest_samples_ls = c()
  bart_sampler_ls  = c()
  
  
  for(mm in 1:M){
    #   outcome_train_ls  = c(outcome_train_ls ,createOutcome(Y[,mm]))
    forest_dataset_train_ls = c(forest_dataset_train_ls, stochtree::createForestDataset(X, basis = NULL, variance_weights = NULL))
    forest_samples_ls = c(forest_samples_ls, stochtree::createForestContainer(num_trees, output_dimension,is_leaf_constant))
    bart_sampler_ls = c(bart_sampler_ls,createForestModel(forest_dataset_train_ls[[mm]], 
                                                          feature_types, num_trees,nrow(X), alpha, beta,
                                                          min_samples_leaf)) 
  }
  
  #Initialization of Covariance Objects:
  th.A0  = matrix(1,M,M)
  eta = matrix(NA, TT, M)
  H = matrix(-3, TT, M)
  d.A0 = diag(M)
  sigma_mat = matrix(NA, M, 1)
  Y.fit.BART = Y*0
  
  
  # Init. for the Horseshoe:
  lambda.A0 = 1
  nu.A0 = 1
  tau.A0 =1
  zeta.A0 =1
  prior.cov = rep(1, M*(M-1)/2)
  
  # Stochastic Volatility:
  sv_priors = list()
  sv_draw = list()
  h_latent = list()
  
  sv_params_mcmc = array(NA, dim = c(num_samples, 4, M))
  sv_params_mat <- matrix(NA, M, 4)
  colnames(sv_params_mat) = c("mu", "phi", "sigma","h")
  
  if(sv){
    for(mm in 1:M){
      sv_draw[[mm]] = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]] = rep(0,TT)
      sv_priors[[mm]] = stochvol::specify_priors(
        mu     = sv_normal(mean =0, sd = 10),
        phi    = sv_beta(shape1 =5 , shape2 = 1.5),
        sigma2 = sv_gamma(shape = 0.5, rate = 10),
        nu     = sv_infinity(),
        rho    = sv_constant(0))
    }
  }else{
    for(mm in 1:M){
      sv_draw[[mm]] = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]] = rep(0,TT)
      sv_priors[[mm]] <- stochvol::specify_priors(
        mu = sv_constant(0),
        phi = sv_constant(1-1e-12),
        sigma2 = sv_constant(1e-12),
        nu = sv_infinity(),
        rho = sv_constant(0)
      )
    }
    
  }
  
  # Container of variance parameter samples
  global_var_samples <- matrix(0, num_samples, M)
  colnames(global_var_samples) = variables
  leaf_scale_samples <- matrix(0, num_samples, M)
  Y.store = array(NA, dim=c(num_samples, TT, M))
  H.store = array(NA, dim=c(num_samples, TT, M))
  # Run MCMC
  
  # Run MCMC
  pb = txtProgressBar(min = 0, max = num_samples, style = 3)
  start = Sys.time()
  
  
  for (i in (num_gfr+1):num_samples) {
    for(mm in 1:M){
      
      if(mm >1){
        eta_mm = eta[,1:(mm -1), drop = FALSE]
        A0_mm = d.A0[mm,1:(mm-1)]
        outcome_train = createOutcome(Y[,mm] -Y.fit.BART[,mm]  - eta_mm%*%A0_mm)
      }else{
        outcome_train = createOutcome(Y[,1] - Y.fit.BART[,1])
      }
      variable_tree_splits = as.integer(c(0,0))
      forest_samples = forest_samples_ls[[mm]] 
      forest_model   = bart_sampler_ls[[mm]]
      forest_dataset_train = forest_dataset_train_ls[[mm]]
      
      forest_model$sample_one_iteration(
        forest_dataset_train, outcome_train, forest_samples, rng,
        feature_types, leaf_model, current_leaf_scale[mm,mm], variable_weights, variable_tree_splits, 
        current_sigma2[mm,mm], cutpoint_grid_size, gfr = FALSE, pre_initialized = FALSE)
      
      leaf_scale_samples[i,mm] <- sample_tau_one_iteration(forest_samples, rng, a_leaf, b_leaf[mm,mm], i-1)
      current_leaf_scale[mm,mm] <- as.matrix(leaf_scale_samples[i,mm])
      global_var_samples[i,mm] <- sample_sigma2_one_iteration(outcome_train, rng, nu, lambda[mm,mm])
      
      current_sigma2[mm,mm] <- global_var_samples[i,mm]
      sigma_mat[mm] = sqrt(global_var_samples[i,mm])
      
      Y.fit.BART[,mm] = forest_samples$predict_raw_single_forest(forest_dataset_train, i-1)
      eta[,mm] = Y[,mm] -  Y.fit.BART[,mm]
      
      if(mm >1){
        norm.mm = as.numeric(exp(-.5*h_latent[[mm]]) * 1/sigma_mat[mm,])
        u_mm    = eta[,1:(mm-1), drop =F]*norm.mm
        eta_mm  = eta[,mm]*norm.mm
        if (mm == 2){
          v0.inv = 1/th.A0[mm,1]
        }else{
          v0.inv = diag(1/th.A0[mm,1:(mm-1)])
        } 
        V.cov = solve(crossprod(u_mm) + v0.inv)
        mu.cov = V.cov %*% crossprod(u_mm, eta_mm)
        mu.cov.draw = mu.cov + t(chol(V.cov)) %*% rnorm(ncol(V.cov)) 
        d.A0[mm,1:(mm-1)] = mu.cov.draw
        
      }
      
    }
    
    
    shocks = eta %*%t(solve(d.A0))
    for (mm in 1:M){
      if(sv == TRUE){
        svdraw_mm =  stochvol::svsample_general_cpp(shocks[,mm]/sigma_mat[mm], 
                                          startpara = sv_draw[[mm]], startlatent = h_latent[[mm]],
                                          priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] = as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        h_latent[[mm]] = svdraw_mm$latent
        sv_params_mat[mm, ] = c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[TT])
        weights = as.numeric(exp(svdraw_mm$latent))
        forest_dataset_train_ls[[mm]] = stochtree::createForestDataset(X, basis = NULL, variance_weights = weights)
        H[,mm] = log(sigma_mat[mm]^2) + svdraw_mm$latent
      }else{
        H[,mm] = log(sigma_mat[mm]^2)
      }
    }
    
    # Updating the Shrinkage priors
    
    hs_draw = get.hs(bdraw = d.A0[lower.tri(d.A0)],lambda.hs = lambda.A0, nu.hs = nu.A0,
                     tau.hs = tau.A0, zeta.hs = zeta.A0)
    lambda.A0 = hs_draw$lambda
    nu.A0  = hs_draw$nu
    tau.A0 = hs_draw$tau
    zeta.A0   = hs_draw$zeta
    prior.cov = hs_draw$psi
    
    th.A0[lower.tri(th.A0)] = prior.cov
    th.A0[th.A0>10] = 10
    th.A0[th.A0<1e-8] = 1e-8
    
    
    H.store[i,,]  = exp(H)
    Y.store[i,,]  = (Y.fit.BART*t(matrix(Ysd,M,TT)))+t(matrix(Ymu,M,TT))
    dimnames(Y.store) = list(paste0("mcmc",1:num_samples),rownames(Y),colnames(Y))
    iter.update = 250
    setTxtProgressBar(pb, i)
    if (i %% iter.update==0){
      end =  Sys.time()
      message(paste0("\n Average time for single draw over last ",iter.update," draws ",
                     round(as.numeric(end-start)/iter.update, digits=4), " seconds, currently at draw ", i))
      start = Sys.time() 
    }
    
  }
  dimnames(Y.store) = list(paste0("mcmc",1:num_samples),rownames(Y),colnames(Y))
  dimnames(H.store) = list(paste0("mcmc",1:num_samples),rownames(Y),colnames(Y))
  return_obj = list("Y.fitted" = Y.store, "H.fitted" = H.store, "Ysd" = Ysd, "Ymu" = Ymu,
                     "sv.mcmc" = sv_params_mcmc, "A0" = d.A0)
  
  return(return_obj)
}


## Sample HS for a regression.

get.hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}






  

    


