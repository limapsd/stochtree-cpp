
#' Run the Multivariate BART Algorithm for Supervised Learning.
#'
#' @param data Outcome to be modeled by the ensemble. This is a TxM matriz.
#' @param lags Number of Lags used in the VAR(p) system.
#' @param bart_prior Which Prior we are going to use in the BART: CGM, DART or Minn.
#' @param sv Stochastic Volatility Flag.
#' @param num_burnin Number of "burn-in" iterations of the MCMC sampler. Default: 100.
#' @param num_mcmc Number of "retained" iterations of the MCMC sampler. Default: 100.
#' @param params 
#' @param num_thin 
#' @param keep_burnin Whether or not "burn-in" samples should be stored. Default FALSE. 
#'
#' @return List of sampling outputs and a wrapper around the sampled forests (which can be used for in-memory prediction on new data).
#' @export
#'
#' @examples
fsv_mbart <- function(data, 
                      lags = 1L,
                      bart_prior = "CGM",
                      sv = "none",
                      num_burnin   =  100L,
                      num_mcmc     = 100L,
                      num_thin     = 1L,
                      keep_burnin = FALSE,
                      params = ls()){
  

###--------------------------------------------------------------------------###
###-------------------- Preprocessing the VAR objects  ----------------------###
###--------------------------------------------------------------------------###

  num_samples <- num_burnin + num_mcmc
  num_saved   <- num_mcmc %/% num_thin
  if(keep_burnin){
    save_set    <- seq(from = num_thin, to = num_samples, by = num_thin)
  }else{
    save_set    <- seq(from = num_thin + num_burnin, to = num_samples, by = num_thin)
  }

  p    <- lags
  M    <- ncol(data) 
  K    <- M*p
  
  Traw <- nrow(data)
  Y_tmp <- as.matrix(data)

  Ymu   <- apply(Y_tmp,  2, mean,na.rm=T)
  Ysd   <- apply(Y_tmp,  2, sd,na.rm=T)
  Y_tmp  <- apply(Y_tmp,  2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})

  # Data checks
  if (any(is.na(Y_tmp))){
    stop("\nNAs in data.\n")
  }
  if (ncol(Y_tmp) < 2) {
    stop("The matrix 'data' should contain at least two variables. \n")
  }
  
  if (is.null(colnames(Y_tmp))) {
    colnames(Y_tmp) <- paste("y", 1:ncol(Y_tmp), sep = "")
    warning(paste("No column names supplied in data, using:",
                  paste(colnames(Y_tmp), collapse = ", "), ", instead.\n"))
  }

  # embed: get lagged values
  X_train <- embed(Y_tmp, dimension = p + 1)[, -(1:M)]
  colnames(X_train) <- paste0(colnames(Y_tmp), ".l", sort(rep(1:p,M)))
  Y_train <- Y_tmp[-c(1:p), ]

  TT <- Traw - p

 if(TT!=nrow(Y_train) | TT!=nrow(X_train)){
   stop("Something went wrong: Tobs != nrow(Y). \n")
 }

  variable_names <- colnames(Y_train)
  feature_types  <- as.integer(rep(0, K))
  
  bart_params <- preprocessmBartParams(params)

  alpha <- bart_params$alpha
  beta  <- bart_params$beta
  
  min_samples_leaf     <- bart_params$min_samples_leaf
  max_depth            <- bart_params$max_depth
  num_trees            <- bart_params$n_trees
  cutpoint_grid_size   <- bart_params$cutpoint_grid_size
  nu_bart              <- bart_params$nu
  lambda_bart          <- bart_params$lambda
  
  a_leaf               <- bart_params$a_leaf
  b_leaf               <- bart_params$b_leaf
  a_forest             <- bart_params$a_forest
  b_forest             <- bart_params$b_forest
  random_seed          <- bart_params$random_seed
  sample_sigma_leaf    <- bart_params$sample_sigma_leaf
  
  alpha_dart           <- bart_params$alpha_dart
  rho_dart             <- bart_params$rho_dart
  a_dart               <- bart_params$a_dart
  b_dart               <- bart_params$b_dart
  
  outcome_model_type <- 0 # numeric
  
###--------------------------------------------------------------------------###
###-------------------------- MLE estimates ---------------------------------###
###--------------------------------------------------------------------------###
  
 XtX <-  t(X_train) %*% X_train
 XtY <-  t(X_train) %*% Y_train
  
  # Initialization of the covariance dependes on the high dimensionality 
  if(TT <= K){
    sigma2hat <- var(Y_train)
  }else{
    beta_ols  <- solve(XtX) %*% XtY
    resid_mle <- Y_train - X_train %*%beta_ols
    sigma2hat <- ( t(resid_mle)%*%(resid_mle) )/TT
  }


 ## Stochtree priors and Initializers:
  quantile_cutoff <- 0.9
  
  if (is.null(lambda_bart)) {
    lambda <- (sigma2hat*qgamma(1-quantile_cutoff,nu_bart))/nu_bart
  }
  
  sigma2_init <- sigma2hat
  if (is.null(b_leaf)) b_leaf <- var(Y_train)/(2*num_trees)
  sigma_leaf_init <- var(Y_train)/(num_trees)
  
  current_leaf_scale <- as.matrix(sigma_leaf_init)
  current_sigma2 <- sigma2_init
  
  
  # Random number generator (std::mt19937)
  if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
  rng <- createRNG(random_seed)

  #Variable Selection Splits
  variable_weight_splits_ls  <- matrix(0, nrow = M, ncol = K) 
  var_count_matrix <- array(NA, dim =  c(num_saved, K, M) )

  
  if (bart_prior == "minn"){
    lag_index <- matrix(0,M,p) # For the minessota prior.  
    for(i in 1:M){
      lag_index[i,] <- seq(i, K, by = M)
    }
  }

  if (bart_prior == "dart"){
    alpha_dart <- 1
    rho_dart <- ncol(X_train)
  } 
  

  # Container of forest samples
  forest_dataset_train_ls <- c()
  forest_container_ls <- c()
  active_forest_ls <- c()
  bart_model_ls     <- c()

  for(mm in 1:M){    
    variable_weight_splits_ls[mm,] <- rep(1/K, K)
    forest_dataset_train_ls <- c(forest_dataset_train_ls, stochtree::createForestDataset(X_train, basis = NULL, variance_weights = NULL) )
    forest_container_ls <- c(forest_container_ls, stochtree::createForestContainer(num_trees, is_leaf_constant =  TRUE))
    bart_model_ls <- c(bart_model_ls, stochtree::createForestModel(forest_dataset_train_ls[[mm]], 
                                                          feature_types, num_trees, nrow(X_train), alpha, beta,
                                                          min_samples_leaf, max_depth))
    active_forest_ls <- c(active_forest_ls, stochtree::createForest(num_trees = num_trees, is_leaf_constant = TRUE))

  }

###--------------------------------------------------------------------------###
###---------------------------- Prior set-up --------------------------------###
###--------------------------------------------------------------------------###


###################################################################################
################################## Covariance Objects #############################
###################################################################################

  # Time varying covariance matrix object that will be used in the forecasting!
  Sig_t <- array(0, dim = c(TT,M,M))
  for(tt in 1:TT) Sig_t[tt,,] = sigma2hat

  # Factor Structure in the Errors
  ut <- Y_train  - X_train%*%beta_ols 
  
  ledermann_bound <- function(m) { as.integer(floor((2*m+1)/2 - sqrt((2*m+1)^2/4 - m^2 + m))) }
  
  Q <- ledermann_bound(M)
  if(Q==1 || Q == 0) Q = 2
  Lambda <-  matrix(0,nrow = M,ncol = Q)
  
  Ft <- prcomp(ut)$x[,1:Q] #Extract the initizaling factors from the residuals. 
  Omega <- matrix(1, TT,Q)
  
  id_Lambda <- matrix(TRUE,M,Q) # Identification, don't think this is necessary...
  theta_Lambda <- Lambda^0
  
  # HS prior for factor loadings (by column)
  lambda_Lambda <- matrix(1,M,Q)
  nu_Lambda <- matrix(1,M,Q)
  zeta_Lambda <- rep(1,Q)
  tau_Lambda  <- rep(1,Q)

  #Initialization of Covariance Objects:
  eta <- matrix(NA, TT, M)
  H <- matrix(-3, TT, M)

  sigma_mat  <- matrix(NA, M, 1)
  Y_fit_BART <- Y_train*0

  # Stochastic Volatility:
  sv_priors <- list()
  sv_draw <- list()
  h_latent <- list()
  
  sv_params_mcmc <- array(NA, dim = c(num_saved, 3, M))
  sv_params_mat  <- matrix(NA, M, 4)
  colnames(sv_params_mat) = c("mu", "phi", "sigma","h")
  
  fsv_priors    <- list()
  fsv_draw      <- list()
  fsv_h_latent  <- list()

  fsv_params_mcmc <- array(NA, dim = c(num_saved, 3, Q))
  fsv_params_mat  <- matrix(NA, Q, 4)
  colnames(fsv_params_mat) <- c("mu", "phi", "sigma","h")

if(sv == "SV"){   
    for(qq in 1:Q){
        fsv_draw[[qq]] <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
        fsv_h_latent[[qq]] <- rep(0,TT)
        fsv_priors[[qq]] <- stochvol::specify_priors(
        mu     = stochvol::sv_normal(mean = 0, sd = 1),
        phi    = stochvol::sv_beta(shape1 = 25, shape2 = 1.5),
        sigma2 = stochvol::sv_gamma(shape = 22, rate = 2.1),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0)) 
    }
    for(mm in 1:M){
      sv_draw[[mm]]   <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]]  <- rep(0,TT)
      sv_priors[[mm]] <- stochvol::specify_priors(
        mu     = stochvol::sv_normal(mean =0, sd = 1),
        phi    = stochvol::sv_beta(shape1 =25 , shape2 = 1.5),
        sigma2 = stochvol::sv_gamma(shape = 22, rate = 2.1),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0))
    }

  }else{
    for(mm in 1:M){
      sv_draw[[mm]]   <- list(mu = 0, phi = 0.99, sigma = 0.01, nu = Inf, rho = 0, beta = NA, latent0 = 0)
      h_latent[[mm]]  <- rep(0,TT)
      sv_priors[[mm]] <- stochvol::specify_priors(
        mu     = stochvol::sv_constant(0),
        phi    = stochvol::sv_constant(1-1e-12),
        sigma2 = stochvol::sv_constant(1e-12),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0)
      )
    } 
  }

  # Container of variance parameter samples
  
  global_var_samples <- matrix(0, num_saved, M)
  colnames(global_var_samples) <- variable_names
  
  leaf_scale_samples <- matrix(0, num_saved, M)

###--------------------------------------------------------------------------###
###---------------------------- Store objects -------------------------------###
###--------------------------------------------------------------------------###
  
  Y_store <- array(NA, dim=c(num_saved, TT, M))
  H_store <- array(NA, dim=c(num_saved, TT, M))
  Sigma_store <- array(NA, dim=c(num_saved, TT, M))
  Hfc_store <- array(NA, dim = c(num_saved, fhorz, M) )
  Loadings_store <-array(NA, dim = c(num_saved, M, Q))
  Y_forecast_store <- array(NA, dim = c(num_saved, fhorz, M) )

  ###--------------------------------------------------------------------------###
###----------------------------- MCMC sampler -------------------------------### 
###----------------------------- for FSV-BAVART -----------------------------### 
###--------------------------------------------------------------------------###
  
  pb = txtProgressBar(min = 0, max = num_samples, style = 3)
  start = Sys.time()
  
  for (i in 1:num_samples) {
    
    is_mcmc <- i > num_burnin
    
    if (is_mcmc) {
      mcmc_counter <- i - (num_burnin)
      if (mcmc_counter %% num_thin == 0) keep_sample <- TRUE
      else keep_sample <- FALSE
    } else {
      if (keep_burnin) keep_sample <- TRUE
      else keep_sample <- FALSE
    }
    
    
    
    
    for(mm in 1:M){  
        # Sampling model Coefficients:
        Y_ <- Y_train - tcrossprod(Ft, Lambda)
        outcome_train <- createOutcome(Y_[,mm] - Y_fit_BART[,mm])

    ###------------------------------ Step 1: -----------------------------------###
    ###---------------------- Sample BART objects--------------------------------###
    ###--------------------------------------------------------------------------###

        #Unpacking the objects for each equation:
        variable_weights      <- variable_weight_splits_ls[mm,]
      
        forest_samples_mean  <- forest_container_ls[[mm]] 
        forest_model_mean    <- bart_model_ls[[mm]]
        active_forest_mean   <- active_forest_ls[[mm]]
        forest_dataset_train <- forest_dataset_train_ls[[mm]]

        forest_model_mean$sample_one_iteration(forest_dataset_train, outcome_train, forest_samples_mean,
                          active_forest_mean, rng, feature_types, outcome_model_type, current_leaf_scale[mm,mm],
                          variable_weights, a_forest, b_forest, current_sigma2[mm,mm], cutpoint_grid_size,
                          keep_forest = keep_sample, gfr = FALSE)

        variable_count_splits <- active_forest_mean$get_forest_split_counts(ncol(X_train))                  
        
        leaf_scale_i <- stochtree::sample_tau_one_iteration(active_forest_mean, rng, a_leaf, b_leaf[mm,mm])
        current_leaf_scale[mm,mm]  <-as.matrix(leaf_scale_i)

        current_sigma2_i <- stochtree::sample_sigma2_one_iteration(outcome_train, rng, nu_bart, lambda_bart[mm,mm])
        current_sigma2[mm,mm] <- current_sigma2_i

        
        
        # After half of the burnin, turn the DART prior on 
        if(i > floor(num_burnin/2)){    
          if(bart_prior == "dart"){
            
            dart_sampler     <- stochtree::sample_dart_splits_one_iteration(variable_count_splits, alpha_dart, rng)
            log_prob_vec     <- dart_sampler$lpv
            variable_weights <- exp(log_prob_vec)
            
            alpha_sampler <- stochtree::sample_alpha_one_iteration(log_prob_vec, a_dart, b_dart, rho_dart, rng)    
            alpha_dart    <- alpha_sampler$alpha
            
          }
          if(bart_prior == "minn"){
            log_probability_spits <- draw_minessota_split(variable_count_splits, mm,lag_index ,diag(sigma2hat),lambda_1, lambda_2, rng)
            variable_weights      <- exp(log_probability_spits)
            
            if(lambda_update == T){
            lambda_1 <- draw_alpha_update(lambda_1, log_probability_spits, 0.5, 1, rho = rho_dart)  
            lambda_2 <- draw_alpha_update(lambda_2, log_probability_spits, 0.5, 1, rho = rho_dart)  
          }
          }
        }  
        variable_weight_splits_ls[mm,] <- variable_weights
        
        
        
        if(i %in% save_set){
          #first time that this appears, so here should start count the index
          index_saved <- which(i == save_set)
          global_var_samples[index_saved, mm] <- current_sigma2_i
          leaf_scale_samples[index_saved, mm] <- leaf_scale_i
          var_count_matrix[index_saved,,mm]   <- variable_count_splits
        }
    
        sigma_mat[mm] <- sqrt(current_sigma2_i)
        
        # i -1 because indexing in cpp starts at 0
        Y_fit_BART[,mm] <- active_forest$predict_raw_single_forest(forest_dataset_train)
        eta[,mm]        <- Y_[,mm] -  Y_fit_BART[,mm]

    }
    ###--------------------- Step 2: Idiossincratic Vol -------------------------###
    ###--------------- Sample stochastic volatility parameters ------------------###
    ###--------------------------------------------------------------------------###

    # norm_mm = as.numeric(exp(-.5*H[,mm])) this look very wrong!
    # eta_mm  = eta[,mm]*norm_mm
    # shocks = eta 
    shocks = Y_train - tcrossprod(Ft, Lambda) - Y_fit_BART
    for (mm in 1:M){
      if(sv == "SV"){
        svdraw_mm =  stochvol::svsample_general_cpp(shocks[,mm]/sigma_mat[mm], 
                                          startpara = sv_draw[[mm]], startlatent = h_latent[[mm]],
                                          priorspec = sv_priors[[mm]])
        sv_draw[[mm]][c("mu", "phi", "sigma")] = as.list(svdraw_mm$para[, c("mu", "phi", "sigma")])
        h_latent[[mm]] = svdraw_mm$latent
        sv_params_mat[mm, ] = c(svdraw_mm$para[, c("mu", "phi", "sigma")], svdraw_mm$latent[TT])
        if(i %in% save_set){
          sv_params_mcmc[index_saved,,mm] = svdraw_mm$para[, c("mu", "phi", "sigma")]
        }

        weights <- as.numeric(exp(svdraw_mm$latent)) # Double Check this later.
        forest_dataset_train_ls[[mm]] <- stochtree::createForestDataset(X_train, basis = NULL, variance_weights = weights)
        H[,mm] <- log(sigma_mat[mm]^2) + svdraw_mm$latent # Sum of the two vols(?)
      }else{
        H[,mm] <- log(sigma_mat[mm]^2)
      }
    }
    H[H<log(1e-6)] <- log(1e-6)
    
    ###------------------------------ Step 3: -----------------------------------###
    ###---------------- Sample the factor loadings and the factors --------------###
    ###--------------------------------------------------------------------------###
    
    ###------------------- Step 3.1: Sample HS hyperparameters ------------------###
    
    resid <- Y_train - Y_fit_BART

    Ft <-  get_factors(resid,S=exp(H),H=exp(Omega),L=Lambda,q=Q,t=TT)
    id.L <- FALSE # whether factor model for the errors should be identified
    if(!id.L) Ft <- apply(Ft, 2, function(x) (x-mean(x))/sd(x)) # normalize factor draw
    Lambda <- get_Lambda(resid, fac = Ft, S = exp(H), pr= theta_Lambda , m=M, q=Q, id.fac=!id_Lambda[1,Q])
    
    ###------------------- Step 3.2: Sample HS hyperparameters ------------------###
    
    # Using a Horseshoe Prior on loadings by colums
    for(qq in 1:Q){
      hs_draw <- get_hs(bdraw=as.numeric(Lambda[id_Lambda[,qq],qq]),
                        lambda.hs=as.numeric(lambda_Lambda[id_Lambda[,qq],qq]),
                        nu.hs=nu_Lambda[id_Lambda[,qq],qq],tau.hs=tau_Lambda[qq],zeta.hs=zeta_Lambda[qq])
      theta_Lambda[id_Lambda[,qq],qq] <- hs_draw$psi # Will be used in the loadings sampling
      
      lambda_Lambda[id_Lambda[,qq],qq] <- hs_draw$lambda
      nu_Lambda[id_Lambda[,qq],qq] <- hs_draw$nu
      tau_Lambda[qq] <- hs_draw$tau
      zeta_Lambda[qq] <- hs_draw$zeta
    }

    ###------------------- Step 3.3: Sample Factor SV ------------------###
    
    if(sv == "SV"){
        for(qq in 1:Q){
            fsvdraw_qq <- stochvol::svsample_general_cpp(Ft[,qq], startpara = fsv_draw[[qq]], 
                                              startlatent = fsv_h_latent[[qq]], priorspec = fsv_priors[[qq]])
            fsv_draw[[qq]][c("mu", "phi", "sigma")] <- as.list(fsvdraw_qq$para[, c("mu", "phi", "sigma")])
            fsv_h_latent[[qq]]   <- fsvdraw_qq$latent
            fsv_params_mat[qq, ] <- c(fsvdraw_qq$para[, c("mu", "phi", "sigma")], fsvdraw_qq$latent[TT])
            if(i %in% save_set){
              fsv_params_mcmc[index_saved,,qq] <- fsvdraw_qq$para[, c("mu", "phi", "sigma")]
            }
            Omega[,qq] <- fsvdraw_qq$latent
      }
    }else{
        for(qq in 1:Q){
            Omega[,qq] <- 0
        }
    }

    for(tt in 1:TT){
      s_t <- Lambda %*%diag( exp(Omega[tt,] ) ) %*% t(Lambda) + diag(exp(H[tt,]))
      Sig_t[tt,,] <- s_t 
    }

    if( i %in% save_set){
      
        ###------------------------------ Step 5: -----------------------------------###
        ###----------------------------- Storage ------------------------------------###
        ###--------------------------------------------------------------------------###

        H_store[index_saved,,]            <- exp(H)
        Y_store[index_saved,,]            <- (Y_fit_BART*t(matrix(Ysd,M,TT)))+t(matrix(Ymu,M,TT))
        Loadings_store[index_saved,,]     <- Lambda
        Omega_store[index_saved,,]        <- Omega
      }

    iter_update <- 250
    setTxtProgressBar(pb, i)
    if (i %% iter_update==0){
      end   <-  Sys.time()
      message(paste0("\n Average time for single draw over last ",iter_update," draws ",
                     round(as.numeric(end-start)/iter_update, digits=4), " seconds, currently at draw ", i))
      start <- Sys.time() 
    }

    }

    model_params <- list(
    "sigma2_init" = sigma2_init, 
    "sigma_leaf_init" = sigma_leaf_init,
    "a_global" = nu_bart,
    "b_global" = lambda_bart, 
    "a_leaf" = a_leaf, 
    "b_leaf" = b_leaf,
    "a_forest" = a_forest, 
    "b_forest" = b_forest,
    "outcome_mean" = Ymu,
    "outcome_scale" = Ysd, 
    "num_covariates" = ncol(X_train), 
    "num_samples" = num_samples, 
    "num_burnin" = num_burnin, 
    "num_mcmc" = num_mcmc, 
    "bart_prior" = bart_prior,
    "sv" = sv 
  )
  result <- list(
    "model_params"           = model_params, 
    "lags"                   = p,
    "y_hat_train"            = Y_store,
    "H_Matrix"               = H_store,
    "L_Matrix"               = Loadings_store,
    "Dt"                     = Omega_store,
    "Yt_train"               = Y_train,
    "Y_raw"                  = data,
    "forest_samples_mean_ls" = forest_container_ls,
    "var_count_matrix"       = var_count_matrix
    )
  
  result[["mean_forests"]] = forest_container_ls
  
  result[["var_count_matrix"]] = var_count_matrix

  return(result)
  }


  

  
  
  


  


  




