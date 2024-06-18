
#' Title
#'
#' @param Y_train 
#' @param X_train 
#' @param Yraw 
#' @param p 
#' @param fhorz 
#' @param Heteroskedacity_Model 
#' @param Variable_Split_Prior 
#' @param Theta_Update 
#' @param cutpoint_grid_size 
#' @param tau_init 
#' @param alpha 
#' @param beta 
#' @param min_samples_leaf 
#' @param output_dimension 
#' @param is_leaf_constant 
#' @param nu 
#' @param leaf_model 
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
varbart =  function(Y_train, X_train, Yraw, p, fhorz, Heteroskedacity_Model = "None",
                    Variable_Split_Prior = "Uniform", Theta_Update = F,   
                    cutpoint_grid_size = 100, tau_init = NULL, alpha = 0.95, 
                    beta = 2.0, min_samples_leaf = 5,  output_dimension = 1,
                    is_leaf_constant = T, nu = 3,leaf_model = 0,
                    lambda = NULL, a_leaf = 3, b_leaf = NULL, q = 0.9,
                    sigma2_init = NULL, num_trees = 250, num_gfr = 0, 
                    num_burnin = 250, num_mcmc = 100, random_seed = -1,
                    keep_burnin = F){
  
  num_samples = num_burnin + num_mcmc 
  
  Ymu = apply(Yraw, 2, mean, na.rm=T)
  Ysd = apply(Yraw, 2, sd, na.rm=T)


  M = ncol(Y_train)
  TT = nrow(Y_train)
  K  = ncol(X_train)
  
  variable_names = colnames(Y_train)
  index_names    = rownames(Y_train)
  feature_types  = rep(0, K)  
  colnames(X_train) = paste(rep(variable_names, p),
                         sort(rep(paste(".t-", sprintf("%02d", 1:p), sep = ""), each = 2)), sep = "" )
  
  #  Calibrate priors for sigma^2 and tau
  
  beta_ols <- solve(crossprod(X_train)) %*% crossprod(X_train,Y_train)
  sigma2hat = crossprod(Y_train - X_train%*%beta_ols)/TT
  
  quantile_cutoff <- 0.9
  if (is.null(lambda)) {
    lambda <- (sigma2hat*qgamma(1-quantile_cutoff,nu))/nu
  }
  if (is.null(sigma2_init)) sigma2_init <- sigma2hat
  if (is.null(b_leaf)) b_leaf <- var(Y_train)/(2*num_trees)
  if (is.null(tau_init)) tau_init <- var(Y_train)/(num_trees)
  current_leaf_scale <- as.matrix(tau_init)
  current_sigma2 <- sigma2_init
  
# Now the time varying varyances, this will be used for the forecasting!
  Sig_t = array(0, dim= c(TT,M,M))
  for(tt in 1:TT) Sig_t[tt,,] = sigma2hat
  

  # Data 
  forest_dataset_train_ls = c()
  
  # Random number generator (std::mt19937)
  if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
  rng <- createRNG(random_seed)
  
  # Sampling data structures
  feature_types <- as.integer(feature_types)
  #outcome_train = createOutcome(Y[,1])
  
  
  # Variable selection weights
    variable_weights <- rep(1/K, K)
    
    #Variable Selection Splits
    variable_count_splits_ls  = matrix(0, nrow = M, ncol = K) 
    var_count_matrix = array(NA, dim =  c(num_samples, K, M) )
    
  
  
  # Container of forest samples
  forest_samples_ls = c()
  bart_sampler_ls  = c()
  
  
  for(mm in 1:M){
    
    variable_count_splits_ls[mm,] =  as.integer(rep(0,K))
    #   outcome_train_ls  = c(outcome_train_ls ,createOutcome(Y[,mm]))
    forest_dataset_train_ls = c(forest_dataset_train_ls, stochtree::createForestDataset(X_train, basis = NULL, variance_weights = NULL) )
    forest_samples_ls = c(forest_samples_ls, stochtree::createForestContainer(num_trees, output_dimension,is_leaf_constant))
    bart_sampler_ls = c(bart_sampler_ls, stochtree::createForestModel(forest_dataset_train_ls[[mm]], 
                                                          feature_types, num_trees, nrow(X_train), alpha, beta,
                                                          min_samples_leaf)) 
  }
  
  #Initialization of Covariance Objects:
  th.A0  = matrix(1,M,M)
  eta = matrix(NA, TT, M)
  H = matrix(-3, TT, M)
  d.A0 = diag(M)
  sigma_mat = matrix(NA, M, 1)
  Y_fit_BART = Y_train*0
  
  #Dirichlet Initialization
  theta = 1
  
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
  

  if(Heteroskedacity_Model == "CSV"){
    for(mm in 1:M){
      sv_draw[[mm]] = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]] = rep(0,TT)
      sv_priors[[mm]] = stochvol::specify_priors(
        mu     = stochvol::sv_normal(mean =0, sd = 10),
        phi    = stochvol::sv_beta(shape1 =5 , shape2 = 1.5),
        sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 10),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0))
    }
  }else{
    for(mm in 1:M){
      sv_draw[[mm]]   = list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]]  = rep(0,TT)
      sv_priors[[mm]] = stochvol::specify_priors(
        mu     = stochvol::sv_constant(0),
        phi    = stochvol::sv_constant(1-1e-12),
        sigma2 = stochvol::sv_constant(1e-12),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0)
      )
    }
    
  }
  
  
  # Container of variance parameter samples
  global_var_samples <- matrix(0, num_samples, M)
  colnames(global_var_samples) = variables
  
  leaf_scale_samples <- matrix(0, num_samples, M)
  Y_store = array(NA, dim=c(num_samples, TT, M))
  H_store = array(NA, dim=c(num_samples, TT, M))
  Y_forecast_store = array(NA, dim = c(num_mcmc, fhorz, M) )

  # Run MCMC
  
  pb = txtProgressBar(min = 0, max = num_samples, style = 3)
  start = Sys.time()
  
  for (i in (num_gfr+1):num_samples) {
    for(mm in 1:M){
      
      if(mm >1){
        eta_mm = eta[,1:(mm -1), drop = FALSE]
        A0_mm = d.A0[mm,1:(mm-1)]
        outcome_train = createOutcome(Y_train[,mm] -Y_fit_BART[,mm]  - eta_mm%*%A0_mm)
      }else{
        outcome_train = createOutcome(Y_train[,1] - Y_fit_BART[,1])
      }
      
      #Unpacking the objects for each equation:
      variable_count_splits = as.integer(variable_count_splits_ls[mm,])
      forest_samples = forest_samples_ls[[mm]] 
      forest_model   = bart_sampler_ls[[mm]]
      forest_dataset_train = forest_dataset_train_ls[[mm]]
      

      variable_count_splits = forest_model$sample_one_iteration(
        forest_dataset_train, outcome_train, forest_samples, rng,
        feature_types, leaf_model, current_leaf_scale[mm,mm], variable_weights, variable_count_splits, 
        current_sigma2[mm,mm], cutpoint_grid_size, gfr = FALSE, pre_initialized = FALSE)

      if(Variable_Split_Prior == "Dirichilet"){
        log_probability_values = draw_probability_splits(variable_count_splits, theta)
        variable_weights = exp(log_probability_values)
        
        if(Theta_Update == TRUE){
          theta = draw_theta_update(theta, log_probability_values, 0.5, 1, rho = length(log_probability_values))  
        }
      }
      var_count_matrix[i,,] = variable_count_splits

      leaf_scale_samples[i,mm] <- sample_tau_one_iteration(forest_samples, rng, a_leaf, b_leaf[mm,mm], i-1)
      current_leaf_scale[mm,mm] <- as.matrix(leaf_scale_samples[i,mm])
      global_var_samples[i,mm] <- sample_sigma2_one_iteration(outcome_train, rng, nu, lambda[mm,mm])
      
      current_sigma2[mm,mm] <- global_var_samples[i,mm]
      sigma_mat[mm] = sqrt(global_var_samples[i,mm])

      # i -1 because indexing in cpp starts at 0
      Y_fit_BART[,mm] = forest_samples$predict_raw_single_forest(forest_dataset_train, i-1)
      eta[,mm] = Y_train[,mm] -  Y_fit_BART[,mm]
      
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
        # d.A0[mm,1:(mm-1)] = mu.cov.draw
        
      }
      
    }
    
    
    shocks = eta %*%t(solve(d.A0))
    for (mm in 1:M){
      if(Heteroskedacity_Model == "CSV"){
        svdraw_mm =  stochvol::svsample_general_cpp(shocks[,mm]/sigma_mat[mm], 
                                          startpara = sv_draw[[mm]], startlatent = h_latent[[mm]],
                                          priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] = as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        h_latent[[mm]] = svdraw_mm$latent
        sv_params_mat[mm, ] = c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[TT])
        
        weights = as.numeric(exp(svdraw_mm$latent)) # Double Check this later.
        forest_dataset_train_ls[[mm]] = stochtree::createForestDataset(X_train, basis = NULL, variance_weights = weights)
        H[,mm] = log(sigma_mat[mm]^2) + svdraw_mm$latent #Is this necessary in the Stochtree Framework? I need to redo this math...
      }else{
        H[,mm] = log(sigma_mat[mm]^2)
      }
    }
    
    for(tt in 1:TT){
      aux_s = exp(H[tt,])
      s_t =  t(d.A0)%*%crossprod(diag(aux_s), d.A0)
      Sig_t[tt,,] = s_t 
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
    
    ################## Creating and Updating the Forecast ####################  
    
    if( i > num_burnin){

      if(fhorz > 0){
        Yfc = matrix(NA,fhorz, M)
        Hfc = matrix(NA, fhorz, M)
        if(p>1){
          X_hat = matrix(NA, nrow =1, ncol = K)
          X_hat[1,] = c(Y_train[TT,], X_train[TT, 1:(M*(p-1))] )
        }else{
          X_hat = Y_train[TT,]
        }
    
        forecast_dataset <- createForestDataset(X_hat)
        if(Heteroskedacity_Model == "CSV"){
          HT = H[TT,] - log(as.numeric(sigma_mat)^2)  
        }else{
          HT = H[TT,]
        }
  
        Sig_T = Sig_t[TT,,]
        forest_predictions = matrix(0, M)
       
        # Run Forecast Loop:
        for(hh in 1:fhorz){
          HT  = log(as.numeric(sigma_mat)^2) + (sv_params_mat[, 1] + sv_params_mat[ , 2] * (HT - sv_params_mat[,1]) + sv_params_mat[ , 3]*rnorm(M))
          Hfc[hh,] = exp(HT)
          
          for(mm in 1:M){
            forest_samples = forest_samples_ls[[mm]] 
            forest_predictions[mm] = forest_samples$predict_raw_single_forest(forecast_dataset, i-1)
          }
          
          if(Heteroskedacity_Model == "CSV"){
            Sig_T = t(d.A0)%*%crossprod(diag(exp(HT)), d.A0)
          }
          Y_tp = as.numeric(forest_predictions) + t(chol(Sig_T)) %*%rnorm(M)
          if(p >1){
            X_hat = matrix(NA, nrow =1, ncol = K)
            X_hat[1,] = c(Y_tp, X_hat[,1: (M*(p-1) )] )
            forecast_dataset <- createForestDataset(X_hat)
          }else{
            X_hat = matrix(NA, nrow =1, ncol = M)
            X_hat[1,] = c(Y_tp)
            forecast_dataset <- createForestDataset(X_hat)
          }
          Yfc[hh,] = Y_tp 
        }
      
        # Y_forecast_store[i,,] = (Yfch*t(matrix(Ysd,M,fhorz)))+t(matrix(Ymu,M,fhorz))
        Y_forecast_store[i - num_burnin,,]   = (Yfc*t(matrix(Ysd,M,fhorz)))+t(matrix(Ymu,M,fhorz))
    }
    
    
  
    
  }

  H_store[i,,]  = exp(H)
  #Y_store[i,,]  = (Y.fit.BART*t(matrix(Ysd,M,TT)))+t(matrix(Ymu,M,TT)) 
  Y_store[i,,]   = Y_fit_BART
  iter_update = 250
  setTxtProgressBar(pb, i)
  if (i %% iter_update==0){
    end =  Sys.time()
    message(paste0("\n Average time for single draw over last ",iter_update," draws ",
                   round(as.numeric(end-start)/iter_update, digits=4), " seconds, currently at draw ", i))
    start = Sys.time() 
  }

  # dimnames(Y_store) = list(paste0("mcmc",1:num_samples),index_names,variable_names)
  # dimnames(H_store) = list(paste0("mcmc",1:num_samples),index_names, variable_names)
  # return_obj = list("Y.fitted" = Y_store, "H.fitted" = H_store, "Ysd" = Ysd, "Ymu" = Ymu,
  #                    "sv.mcmc" = sv_params_mcmc, "A0" = d.A0)
  
  
}

  dimnames(Y_store) = list(paste0("mcmc", 1:num_samples), index_names, variable_names)
  dimnames(H_store) = list(paste0("mcmc",1:num_samples),index_names, variable_names)

  model_params <- list(
      "lags" = p,
      "forecast_horizon" = fhorz,
      "Variable_Split_Prior" = Variable_Split_Prior,
      "Heteroskedacity_Model" = Heteroskedacity_Model,
      "sigma2_init" = sigma2_init, 
      "nu" = nu,
      "lambda" = lambda, 
      "tau_init" = tau_init,
      "a" = a_leaf, 
      "b" = b_leaf,
      "outcome_mean" = Ymu,
      "outcome_scale" = Ysd, 
      "output_dimension" = output_dimension,
      "is_leaf_constant" = is_leaf_constant,
      "num_endogenous_variables" = M,
      "num_covariates" = K,  
      "num_samples" = num_samples, 
      "num_gfr" = num_gfr, 
      "num_burnin" = num_burnin, 
      "num_mcmc" = num_mcmc, 
      "sv_params_mcmc" = sv_params_mcmc
  )
  result <- list(
      "forests" = forest_samples_ls, 
      "model_params" = model_params, 
      "y_hat_insample" = Y_store,
      "y_forecast"     = Y_forecast_store,
      "H_matrix"       = H_store,
      "variable_count_splits" = var_count_matrix,
      "A0_matrix" = d.A0
  )

 return(result)
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






  

    


