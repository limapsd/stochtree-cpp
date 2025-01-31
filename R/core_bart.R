
#' Run the BART algorithm for supervised learning with the basic CGM priors and Sparse option.
#'
#' @param y_train Outcome to be modeled by the ensemble.
#' @param X_train Covariates used to split trees in the ensemble. May be provided either as a dataframe or a matrix. 
#' Matrix covariates will be assumed to be all numeric. Covariates passed as a dataframe will be 
#' preprocessed based on the variable types (e.g. categorical columns stored as unordered factors will be one-hot encoded, 
#' categorical columns stored as ordered factors will passed as integers to the core algorithm, along with the metadata 
#' that the column is ordered categorical).
#' @param X_test (Optional) Test set of covariates used to define "out of sample" evaluation data. 
#' May be provided either as a dataframe or a matrix, but the format of `X_test` must be consistent with 
#' that of `X_train`.
#' @param alpha Prior probability of splitting for a tree of depth 0 in the mean model. Tree split prior combines `alpha_mean` and `beta_mean` via `alpha_mean*(1+node_depth)^-beta_mean`. Default: `0.95`.
#' @param beta Exponent that decreases split probabilities for nodes of depth > 0 in the mean model. Tree split prior combines `alpha_mean` and `beta_mean` via `alpha_mean*(1+node_depth)^-beta_mean`. Default: `2`.
#' @param min_samples_leaf Minimum allowable size of a leaf, in terms of training samples, in the mean model. Default: `5`.
#' @param max_depth Maximum depth of any tree in the ensemble in the mean model. Default: `10`. Can be overridden with ``-1`` which does not enforce any depth limits on trees.
#' @param n_trees Number of trees in the ensemble for the conditional mean model. Default: `200`. If `n_trees = 0`, the conditional mean will not be modeled using a forest, and the function will only proceed if `num_trees_variance > 0`.
#' @param cutpoint_grid_size Maximum size of the "grid" of potential cutpoints to consider. Default: `100`.
#' @param outcome_train 
#' @param nu 
#' @param lambda 
#' @param a_leaf Shape parameter in the `IG(a_leaf, b_leaf)` leaf node parameter variance model. Default: `3`.
#' @param b_leaf Scale parameter in the `IG(a_leaf, b_leaf)` leaf node parameter variance model. Calibrated internally as `0.5/n_trees` if not set here.
#' @param a_forest 
#' @param b_forest 
#' @param num_burnin Number of "burn-in" iterations of the MCMC sampler. Default: 100.
#' @param num_mcmc Number of "retained" iterations of the MCMC sampler. Default: 100.
#' @param sample_sigma_global 
#' @param sample_sigma_leaf 
#' @param sparse
#' @param alpha_dart 
#' @param rho_dart 
#' @param a_dart 
#' @param b_dart 
#'
#' @return List of sampling outputs and a wrapper around the sampled forests
#' @export
#'
#' @examples
core_bart <-function(y_train, X_train,
                     X_test = NULL,
                     alpha = 0.95,
                     beta = 2.0,  
                     min_samples_leaf = 5,
                     max_depth = 10,
                     n_trees = 200,
                     cutpoint_grid_size = 100,
                     outcome_train = 1,
                     nu = 4,
                     lambda = 0.5,
                     a_leaf = 2,
                     b_leaf = NULL,
                     a_forest = 1,
                     b_forest = 1,
                     num_burnin = 100,
                     num_mcmc = 100,
                     random_seed = -1,
                     sample_sigma_global = TRUE,
                     sample_sigma_leaf   = TRUE,
                     sparse = FALSE,
                     alpha_dart = NULL, 
                     rho_dart = NULL,
                     a_dart = 1,
                     b_dart = 0.5,
                     keep_burnin = FALSE){
  
  
  num_samples <- num_burnin + num_mcmc
  
  # Preprocess covariates
  if ((!is.data.frame(X_train)) && (!is.matrix(X_train))) {
    stop("X_train must be a matrix or dataframe")
  }
  
  if (!is.null(X_test)){
    if ((!is.data.frame(X_test)) && (!is.matrix(X_test))) {
      stop("X_test must be a matrix or dataframe")
    }
  }
  
  train_cov_preprocess_list <- preprocessTrainData(X_train)
  X_train_metadata <- train_cov_preprocess_list$metadata
  X_train <- train_cov_preprocess_list$data
  original_var_indices <- X_train_metadata$original_var_indices
  feature_types <- X_train_metadata$feature_types
  if (!is.null(X_test)) X_test <- preprocessPredictionData(X_test, X_train_metadata)
  
  # Data preprocessing:
  y_bar_train <- mean(y_train)
  y_std_train <- sd(y_train)
  resid_train <- (y_train-y_bar_train)/y_std_train
  
  
  # Convert y_train to numeric vector if not already converted
  if (!is.null(dim(y_train))) {
    y_train <- as.matrix(y_train)
  }
  
  # Determine whether a test set is provided
  has_test = !is.null(X_test)
  
  # Stochtree initializers:
  sigma_leaf_init <- var(resid_train)/(n_trees)
  current_leaf_scale <- as.matrix(sigma_leaf_init)
  #   feature_types <- as.integer(rep(0, ncol(X_train))) # 0 = numeric
  variable_weights_mean <- rep(1/ncol(X_train), ncol(X_train))
  
  # Update variable weights
  variable_weights_adj <- 1/sapply(original_var_indices, function(x) sum(original_var_indices == x))
  variable_weights_mean <- variable_weights_mean[original_var_indices]*variable_weights_adj
  
  
  outcome_model_type <- 0 # numeric
  outcome_train <- stochtree::createOutcome(resid_train)
  
  if (sparse){
    alpha_dart <- 1
    rho_dart <- ncol(X_train)
  }
  
  # Calibrate priors for sigma^2 and tau
  if (is.null(b_leaf)) b_leaf <- var(resid_train)/(2*n_trees)
  sigma2_init <- var(resid_train)  
  current_sigma2 <- sigma2_init  
  # Store objects:
  sigma2_samples     <- c(sigma2_init, rep(0, num_samples))
  leaf_scale_samples <- c(current_leaf_scale, rep(0, num_samples))
  var_count_matrix   <- matrix(NA, nrow = num_samples, ncol = ncol(X_train)) 
  
  # Random number generator (std::mt19937)
  if (is.null(random_seed)) random_seed = sample(1:10000,1,F)
  rng <- stochtree::createRNG(random_seed)
  
  # Creating Stochtree objects:
  forest_dataset_mean <- stochtree::createForestDataset(X_train) # store dataset
  if(has_test) forest_dataset_test <- stochtree::createForestDataset(X_test) 
  
  # Sampling data structures:
  feature_types <- as.integer(feature_types)
  
  forest_model_mean   <- stochtree::createForestModel(forest_dataset_mean, feature_types, 
                                                      n_trees, nrow(X_train), alpha, beta,
                                                      min_samples_leaf, max_depth)
  
  forest_samples_mean <-stochtree::createForestContainer(num_trees = n_trees, is_leaf_constant = T)
  active_forest_mean  <-stochtree::createForest(num_trees = n_trees, is_leaf_constant = T)
  
  
  pb = txtProgressBar(min = 0, max = num_samples, style = 3)
  start = Sys.time()
  
  # Start Running the MCMC
  
  for(i in 1:num_samples){
    
    is_mcmc <- i > num_burnin
    
    if (is_mcmc) {
      mcmc_counter <- i - (num_burnin)
      if (mcmc_counter %% 1 == 0) keep_sample <- TRUE
      else keep_sample <- FALSE
    } else {
      if (keep_burnin) keep_sample <- TRUE
      else keep_sample <- FALSE
    }
    
    
    forest_model_mean$sample_one_iteration(
      forest_dataset_mean, outcome_train, forest_samples_mean, active_forest_mean, rng, feature_types,
      outcome_model_type, current_leaf_scale, variable_weights_mean, a_forest, b_forest, current_sigma2,
      cutpoint_grid_size, keep_forest = keep_sample, gfr = FALSE)
    
    variable_count_splits <- active_forest_mean$get_forest_split_counts(ncol(X_train))
    var_count_matrix[i,]  <- variable_count_splits
    
    # Sample global variance parameter
    if (sample_sigma_global){
      current_sigma2  <- stochtree::sample_sigma2_one_iteration(outcome_train, 
                                                                forest_dataset_mean,
                                                                rng, nu, lambda)
      sigma2_samples[i+1]<- current_sigma2
    }
    
    # Sample leaf node variance parameter and update `leaf_prior_scale`
    if (sample_sigma_leaf){
      leaf_scale_samples[i+1] <- stochtree::sample_tau_one_iteration(active_forest_mean,
                                                                     rng, a_leaf, b_leaf)
      current_leaf_scale[1,1] <- leaf_scale_samples[i+1]
    }
    
    if (sparse){
      if(i > floor(num_burnin/2) ){
    
        dart_sampler      <- sample_dart_splits_one_iteration(variable_count_splits, alpha_dart, rng)
        log_prob_vec      <- dart_sampler$lpv
        variable_weights_mean <- exp(log_prob_vec)
        
        alpha_sampler <- sample_alpha_one_iteration(log_prob_vec, a_dart, b_dart, rho_dart, rng)
        alpha_dart    <- alpha_sampler$alpha
      }
    }
    
    iter_update = 250
    setTxtProgressBar(pb, i)
    if (i %% iter_update==0){
      end =  Sys.time()
      message(paste0("\n Average time for single draw over last ",iter_update," draws ",
                     round(as.numeric(end-start)/iter_update, digits=4), " seconds, currently at draw ", i))
      start = Sys.time() 
    }
  }
  
  # Mean forest predictions
  y_hat_train <- forest_samples_mean$predict(forest_dataset_mean)*y_std_train + y_bar_train
  if (has_test) y_hat_test <- forest_samples_mean$predict(forest_dataset_test)*y_std_train + y_bar_train
  
  
  
  # Prepping return objects:
  if (sample_sigma_global) {
    sigma2_samples <- sigma2_samples[(num_burnin+1):length(sigma2_samples)]
  }
  if (sample_sigma_leaf) {
    leaf_scale_samples <- leaf_scale_samples[(num_burnin+1):length(leaf_scale_samples)]
  }
  
  var_count_matrix <- var_count_matrix[(num_burnin+1):num_samples,]
  
  model_params <- list(
    "sigma2_init" = sigma2_init, 
    "sigma_leaf_init" = sigma_leaf_init,
    "a_global" = nu,
    "b_global" = lambda, 
    "a_leaf" = a_leaf, 
    "b_leaf" = b_leaf,
    "a_forest" = a_forest, 
    "b_forest" = b_forest,
    "outcome_mean" = y_bar_train,
    "outcome_scale" = y_std_train, 
    "num_covariates" = ncol(X_train), 
    "num_samples" = num_samples, 
    "num_burnin" = num_burnin, 
    "num_mcmc" = num_mcmc, 
    "sample_sigma_global" = sample_sigma_global,
    "sample_sigma_leaf" = sample_sigma_leaf,
    "sparse" = sparse
  )
  result <- list(
    "model_params" = model_params, 
    "train_set_metadata" = X_train_metadata
  )
  
  result[["mean_forests"]] = forest_samples_mean
  result[["y_hat_train"]] = y_hat_train
  result[["var_count_matrix"]] = var_count_matrix
  if (has_test) result[["y_hat_test"]] = y_hat_test
  
  rm(forest_model_mean)
  rm(forest_dataset_mean)
  if (has_test) rm(forest_dataset_test)
  rm(outcome_train)
  rm(rng)
  
  return(result)
  }
  