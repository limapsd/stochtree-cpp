#------------- function to run BART + StochVol -------------#
#' Run the BART algorithm for supervised learning with the basic CGM priors and Sparse option + SV
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
#' @param num_burnin Number of "burn-in" iterations of the MCMC sampler. Default: 1000L.
#' @param num_mcmc Number of "retained" iterations of the MCMC sampler. Default: 1000L.
#' @param num_thin 
#' @param general_params 
#' @param mean_forest_params 
#' @param mean_prior 
#' @param SV 
#' @param sv_priors_list 
#'
#' @return
#' @export
#'
#' @examples
bart_cgm <- function(
    y_train,
    X_train,
    X_test = NULL,
    num_burnin = 1000L,
    num_mcmc = 1000L,
    num_thin = 1L,
    general_params = list(),
    mean_forest_params = list(),
    mean_prior = "cgm", #c("cgm", "dart")
    SV = FALSE,
    sv_priors_list = list("priormu" = c(0, 1), "priorphi" = c(10, 3), "priorsigma" = c(0.5,0.5))
) {
  
  # --- General params ---
  general_params_default <- list(
    cutpoint_grid_size = 100,
    standardize = TRUE,
    sample_sigma2_global = TRUE,
    sigma2_global_init = NULL,
    sigma2_global_shape = 1,
    sigma2_global_scale = 1,
    variable_weights = NULL,
    random_seed = -1,
    keep_burnin = FALSE,
    keep_gfr = FALSE,
    keep_every = 1,
    num_chains = 1,
    verbose = FALSE,
    num_threads = -1
  )
  general_params_updated <- preprocessParams(general_params_default, general_params)
  
  # --- Mean-forest (added τ knobs) ---
  mean_forest_params_default <- list(
    num_trees = 200,
    alpha = 0.95,
    beta = 2.0,
    min_samples_leaf = 5,
    max_depth = 10,
    fixed_sigma2_leaf  = FALSE,
    sample_sigma2_leaf = TRUE,     # τ sampling (leaf prior variance)
    sigma2_leaf_init = NULL,        # init τ
    sigma2_leaf_shape = 3,       
    sigma2_leaf_scale = NULL,    
    keep_vars = NULL,
    drop_vars = NULL,
    num_features_subsample = NULL,
    lags = 4,
    a_dart = 0.5,
    b_dart = 1.
  )
  mean_forest_params_updated <- preprocessParams(mean_forest_params_default, mean_forest_params)
  
  # --- Unpack ---
  cutpoint_grid_size   <- general_params_updated$cutpoint_grid_size
  sample_sigma2_global <- general_params_updated$sample_sigma2_global
  sigma2_init          <- general_params_updated$sigma2_global_init
  a_global             <- general_params_updated$sigma2_global_shape
  b_global             <- general_params_updated$sigma2_global_scale
  variable_weights     <- general_params_updated$variable_weights
  random_seed          <- general_params_updated$random_seed
  keep_burnin          <- general_params_updated$keep_burnin
  keep_gfr             <- general_params_updated$keep_gfr
  keep_every           <- general_params_updated$keep_every
  num_chains           <- general_params_updated$num_chains
  verbose              <- general_params_updated$verbose
  num_threads          <- general_params_updated$num_threads
  
  num_trees_mean       <- mean_forest_params_updated$num_trees
  alpha_mean           <- mean_forest_params_updated$alpha
  beta_mean            <- mean_forest_params_updated$beta
  min_samples_leaf_mean<- mean_forest_params_updated$min_samples_leaf
  max_depth_mean       <- mean_forest_params_updated$max_depth
  sample_sigma2_leaf   <- mean_forest_params_updated$sample_sigma2_leaf
  sigma2_leaf_init     <- mean_forest_params_updated$sigma2_leaf_init
  a_leaf               <- mean_forest_params_updated$sigma2_leaf_shape
  b_leaf               <- mean_forest_params_updated$sigma2_leaf_scale
  keep_vars_mean       <- mean_forest_params_updated$keep_vars
  drop_vars_mean       <- mean_forest_params_updated$drop_vars
  lags                 <- mean_forest_params_updated$lags
  fixed_sigma2_leaf    <- mean_forest_params_updated$fixed_sigma2_leaf
  a_dart               <- mean_forest_params_updated$a_dart
  b_dart               <- mean_forest_params_updated$b_dart
  
  
  # --- Schedule ---
  num_samples <- num_burnin + num_mcmc
  num_saved   <- num_mcmc %/% num_thin
  if (keep_burnin) {
    save_set <- seq(from = num_thin, to = num_samples, by = num_thin)
  } else {
    save_set <- seq(from = num_burnin + num_thin, to = num_samples, by = num_thin)
  }
  if (is.null(variable_weights)) {
    variable_weights_mean <- rep(1 / ncol(X_train), ncol(X_train))
  }
  
  # --- Data preprocessing ---
  train_cov_preprocess_list <- preprocessTrainData(X_train)
  X_train_metadata <- train_cov_preprocess_list$metadata
  X_train <- train_cov_preprocess_list$data
  original_var_indices <- X_train_metadata$original_var_indices
  feature_types <- X_train_metadata$feature_types
  if (!is.null(X_test)) X_test <- preprocessPredictionData(X_test, X_train_metadata)
  
  variable_weights_adj <- 1 / sapply(original_var_indices, function(x) sum(original_var_indices == x))
  variable_weights_mean <- variable_weights_mean[original_var_indices] * variable_weights_adj
  
  if ((!is.null(X_test)) && (ncol(X_test) != ncol(X_train))) {
    stop("X_train and X_test must have the same number of columns")
  }
  if (!is.null(dim(y_train))) y_train <- as.matrix(y_train)
  has_test <- !is.null(X_test)
  
  y_bar_train  <- mean(y_train)
  y_std_train  <- sd(y_train)
  resid_train  <- (y_train - y_bar_train) / y_std_train
  
  # --- AR(p) init for global sigma2 ---
  data_arp <- embed(resid_train, lags + 1)
  yt_arp   <- data_arp[, 1]
  Xt_arp   <- data_arp[, -1]
  fit_arp  <- lm(yt_arp ~ Xt_arp)
  sigma2_init <- (summary(fit_arp)$sigma)^2
  current_sigma2 <- sigma2_init
  
  
  if (is.null(b_leaf)) b_leaf <- var(resid_train)/ (2 * num_trees_mean)
  
  # init τ at its prior mean tau0 unless user provided one
  if (is.null(sigma2_leaf_init)) {
    sigma2_leaf_init <- as.matrix( 2 * var(resid_train) / (num_trees_mean) )
  }
  current_leaf_scale <- sigma2_leaf_init
  
  if (fixed_sigma2_leaf) {
    k_hyper <- 2  
    sample_sigma2_leaf <- FALSE
    
    ymin <- min(resid_train)
    ymax <- max(resid_train)
    
    node_scale <- (ymax - ymin) / 2       
    sigma_mu   <- node_scale / (k_hyper * sqrt(num_trees_mean))
    current_leaf_scale <- sigma_mu^2      
  }
  
  if(mean_prior == "dart"){
    alpha_dart <-1.
    rho_dart <- ncol(X_train)
  }
  
  # --- SV setup (idiosyncratic) ---
  if (SV) {
    Tobs    <- length(resid_train)
    sv_params_mcmc <- matrix(
      NA_real_, nrow = num_saved, ncol = 3L,
      dimnames = list(paste0("mcmc_", seq_len(num_saved)), c("mu","phi","sigma"))
    )
    
    sv_params <- list(
      mu      = 0 + rnorm(1, sd = 0.1), # Since resid_train is standardized
      phi     = 0.8 + pmin(rnorm(1, sd = 0.06), 0.095),
      sigma   = 0.1 + rgamma(1, shape = 1, rate = 10),
      nu      = Inf,
      rho     = 0,
      beta    = NA,
      latent0 = 0
    )
    sv_params$latent0 <- sv_params$mu + rnorm(1)
    h_t <- sv_params$mu + rnorm(Tobs)
    
    if(is.null(sv_priors_list)){
    sv_priors <- stochvol::specify_priors(
      mu     = stochvol::sv_normal(mean = 0, sd = 1), #Since data is standardized, the prior on the uncondicional mean should reflect that.
      phi    = stochvol::sv_beta(shape1 = 25, shape2 = 1.5),
      sigma2 = stochvol::sv_inverse_gamma(shape = 22, rate = 2.1),
      nu     = stochvol::sv_infinity(),
      rho    = stochvol::sv_constant(0)
    )
    }else{
      prior_mu<-sv_priors_list$priormu
      prior_phi<-sv_priors_list$priorphi
      prior_sigma<-sv_priors_list$priorsigma
      
      sv_priors <- stochvol::specify_priors(
        mu     = stochvol::sv_normal(mean = prior_mu[1], sd = prior_mu[2]), #Since data is standardized, the prior on the uncondicional mean should reflect that.
        phi    = stochvol::sv_beta(shape1 = prior_phi[1], shape2 = prior_phi[2]),
        sigma2 = stochvol::sv_gamma(shape = prior_sigma[1], rate = prior_sigma[2]),
        nu     = stochvol::sv_infinity(),
        rho    = stochvol::sv_constant(0)
      )
      
    }
    current_sigma2 <- 1.0
  }
  
  # --- Model objects ---
  leaf_model_mean_forest <- 0
  leaf_dimension <- 1
  is_leaf_constant <- TRUE
  leaf_regression <- FALSE
  init_values_mean_forest <- 0.
  
  forest_dataset_train <- stochtree::createForestDataset(X_train, basis = NULL, variance_weights = rep(1, nrow(X_train)))
  if (has_test) forest_dataset_test <- stochtree::createForestDataset(X_test)
  requires_basis <- FALSE
  outcome_train <- stochtree::createOutcome(resid_train)
  
  if (is.null(random_seed)) random_seed <- sample(1:10000, 1, FALSE)
  rng <- stochtree::createCppRNG(random_seed)
  feature_types <- as.integer(feature_types)
  global_model_config <- stochtree::createGlobalModelConfig(global_error_variance = current_sigma2)
  
  forest_model_config_mean <- stochtree::createForestModelConfig(
    feature_types = feature_types,
    num_trees = num_trees_mean,
    num_features = ncol(X_train),
    num_observations = nrow(X_train),
    variable_weights = variable_weights_mean,
    leaf_dimension = leaf_dimension,
    alpha = alpha_mean,
    beta = beta_mean,
    min_samples_leaf = min_samples_leaf_mean,
    max_depth = max_depth_mean,
    leaf_model_type = leaf_model_mean_forest,
    leaf_model_scale = current_leaf_scale,   
    cutpoint_grid_size = cutpoint_grid_size,
    num_features_subsample = NULL
  )
  forest_model_mean <- stochtree::createForestModel(
    forest_dataset_train, forest_model_config_mean, global_model_config
  )
  forest_samples_mean <- stochtree::createForestSamples(num_trees_mean, leaf_dimension, is_leaf_constant, FALSE)
  active_forest_mean  <- stochtree::createForest(num_trees_mean, leaf_dimension, is_leaf_constant, FALSE)
  
  # --- Storage ---
  sigma2_leaf_samples   <- rep(0, num_saved)
  sigma2_samples   <- rep(0, num_saved)
  var_count_matrix <- matrix(NA, nrow = num_saved, ncol = ncol(X_train))
  var_prob_matrix  <- matrix(NA, nrow = num_saved, ncol = ncol(X_train))
  y_hat_store      <- matrix(NA, nrow = num_saved, ncol = nrow(X_train))
  if (SV) sigt_store <- matrix(NA, nrow = num_saved, ncol = nrow(X_train))
  
  pb <- txtProgressBar(min = 0, max = num_samples, style = 3)
  start <- Sys.time()
  
  active_forest_mean$prepare_for_sampler(
    forest_dataset_train, outcome_train, forest_model_mean,
    leaf_model_mean_forest, init_values_mean_forest
  )
  
  for (i in 1:num_samples) {
    is_mcmc <- i > num_burnin
    if (is_mcmc) {
      mcmc_counter <- i - num_burnin
      keep_sample <- if (has_test) FALSE else (mcmc_counter %% num_thin == 0)
    } else {
      keep_sample <- if (has_test || !keep_burnin) FALSE else TRUE
    }
    
    forest_model_mean$sample_one_iteration(
      forest_dataset = forest_dataset_train,
      residual = outcome_train,
      forest_samples = forest_samples_mean,
      active_forest = active_forest_mean,
      rng = rng,
      forest_model_config = forest_model_config_mean,
      global_model_config = global_model_config,
      num_threads = num_threads,
      keep_forest = keep_sample,
      gfr = FALSE
    )
    
    mu_fit <- active_forest_mean$predict(forest_dataset_train)
    
    if (!SV) {
      current_sigma2 <- stochtree::sampleGlobalErrorVarianceOneIteration(
        outcome_train, forest_dataset_train, rng, a_global, b_global)
      global_model_config$update_global_error_variance(current_sigma2)
    } else {
        shocks <- (resid_train - mu_fit)
        svdraw <- stochvol::svsample_general_cpp(
          shocks, startpara = sv_params, startlatent = h_t, priorspec = sv_priors)
        
        sv_params[c("mu","phi","sigma")] <- as.list(svdraw$para[, c("mu","phi","sigma")])
        h_raw <- as.numeric(svdraw$latent)
        mu_h  <- as.numeric(svdraw$para[, "mu"])
        current_sigma2 <- exp(mu_h)
        u_t <- h_raw - mu_h
        u_t <- pmax(u_t, log(1e-6))
        wts <- exp(u_t)
        h_t <- h_raw
      
      forest_dataset_train$update_variance_weights(variance_weights = wts)
      global_model_config$update_global_error_variance(current_sigma2)
    }
    
    # ------- τ (leaf prior variance) update with cap & schedule -------
    if (sample_sigma2_leaf) {
        tau_draw <- sampleLeafVarianceOneIteration(
          active_forest_mean, rng, a_leaf, b_leaf
        )
        current_leaf_scale <- as.matrix(tau_draw)
        forest_model_config_mean$update_leaf_model_scale(current_leaf_scale)
      }
    
    variable_count_splits <- active_forest_mean$get_forest_split_counts(ncol(X_train))

    if (i > floor(num_burnin / 2)) {
        if(mean_prior == "dart"){
            draw_dart <- sample_dart_splits_one_iteration(variable_count_splits, alpha_dart, rng)
            log_prob_vector <- draw_dart$lpv
            variable_prob_splits <- exp(log_prob_vector)
            
            
            alpha_sampler <- sample_alpha_one_iteration(log_prob_vector, a_dart, b_dart, rho_dart, rng)
            alpha_dart    <- alpha_sampler$alpha
            
            forest_model_config_mean$update_variable_weights(exp(log_prob_vector))
        }
      
    }
    
    
    if (i %in% save_set) {
      idx <- which(save_set == i)
      y_hat_store[idx, ] <- mu_fit * y_std_train + y_bar_train
      sigma2_samples[idx] <- current_sigma2 * (y_std_train ^ 2)
      sigma2_leaf_samples[idx] <- current_leaf_scale
      var_count_matrix[idx, ] <- variable_count_splits
      if(mean_prior =="dart") var_prob_matrix[idx,] <- variable_prob_splits
      if (SV) {
        sigt_store[idx, ] <- exp(h_t/2) * (y_std_train)
        sv_params_mcmc[idx, ] <- unlist(sv_params[c("mu","phi","sigma")])
      }
    }
    
    iter_update <- 250
    setTxtProgressBar(pb, i)
    if (i %% iter_update == 0) {
      end <- Sys.time()
      message(paste0(
        "\n Average time for single draw over last ", iter_update, " draws ",
        round(as.numeric(end - start)/iter_update, 4),
        " seconds, currently at draw ", i
      ))
      start <- Sys.time()
    }
  } # end MCMC loop
  close(pb)
  
  model_params <- list(
    "sigma2_init" = sigma2_init,
    "a_global" = a_global,
    "b_global" = b_global,
    "a_leaf" = a_leaf,
    "b_leaf" = b_leaf,
    "outcome_mean" = y_bar_train,
    "outcome_scale" = y_std_train,
    "leaf_dimension" = leaf_dimension,
    "is_leaf_constant" = is_leaf_constant,
    "leaf_regression" = leaf_regression,
    "requires_basis" = requires_basis,
    "num_covariates" = ncol(X_train),
    "num_burnin" = num_burnin,
    "num_mcmc" = num_mcmc,
    "sample_sigma2_global" = sample_sigma2_global
  )
  
  result <- list("model_params" = model_params, "train_set_metadata" = X_train_metadata)
  result[["y_hat_train"]] <- y_hat_store
  if (has_test) result[["y_hat_test"]] <- y_hat_test_store
  if (SV) {
    result[["sigmat_store"]] <- sigt_store
    result[["sv_params_mcmc"]] <- sv_params_mcmc
  } else {
    result[["sigma2_draws"]] <- sigma2_samples
  }
  result[["sigma2_leaf_draws"]] <- sigma2_leaf_samples
  result[["var_count_matrix"]] <- var_count_matrix
  result[["var_prob_matrix"]] <- var_prob_matrix
  class(result) <- "bartmodel"
  
  rm(forest_model_mean); rm(forest_dataset_train); rm(outcome_train); rm(rng)
  if (has_test) rm(forest_dataset_test)
  
  return(result)
}
