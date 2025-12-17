#' Title
#'
#' @param data 
#' @param Y_test 
#' @param n_ahead 
#' @param lags 
#' @param mean_prior 
#' @param variance_prior 
#' @param cov_shrinkage 
#' @param num_factors 
#' @param num_warmstart 
#' @param num_burnin 
#' @param num_mcmc 
#' @param num_thin 
#' @param general_params 
#' @param mean_forest_params 
#'
#' @return
#' @export
#'
#' @examples
varbart <- function(
    data,
    Y_test = NULL,
    n_ahead = NULL,
    lags = 1L,
    mean_prior = "cgm", # c("cgm", "sparse","minnesota")
    variance_prior = "const", # c("const", "csv", "fsv")
    cov_shrinkage = "hs", # c("hs","ng") for the factor loadings
    num_factors = NULL,
    num_warmstart = 0L,
    num_burnin = 100L,
    num_mcmc = 100L,
    num_thin = 1L,
    general_params = list(),
    mean_forest_params = list()
) {
  
    ###--------------------------------------------------------------------------###
    ###-------------------- Preprocessing the VAR objects  ----------------------###
    ###--------------------------------------------------------------------------###
    
    p <- lags
    M <-ncol(data)
    
    if( M == 1 ){
      stop("For a VAR use a T x M matrix with M >= 2.")
    }

    if( mean_prior != "cgm" && mean_prior != "sparse" && mean_prior != "minnesota") {
      stop("mean_prior must be one of 'cgm', 'sparse' or 'minnesota'.")
    }

    if( cov_shrinkage != "hs" && cov_shrinkage != "ng") {
      stop("cov_shrinkage must be one of 'hs' or 'ng'.")
    }

    if (!variance_prior %in% c("const","csv","fsv")) {
      stop("variance_prior must be one of 'const', 'csv' or 'fsv'.")
    }

    K <- M * p

    Traw <- nrow(data)
    Y_tmp <- as.matrix(data)
    
    Ymu <- apply(Y_tmp, 2, mean, na.rm = T)
    Ysd <- apply(Y_tmp, 2, sd, na.rm = T)
    
    Y_tmp <- apply(Y_tmp, 2, function(x) {
      (x - mean(x, na.rm = T)) / sd(x, na.rm = T)
    })
    
    # Data checks
    if (any(is.na(Y_tmp))) {
      stop("\nNAs in data.\n")
    }
       
    if (is.null(colnames(Y_tmp))) {
      colnames(Y_tmp) <- paste("y", 1:ncol(Y_tmp), sep = "")
      warning(paste(
        "No column names supplied in data, using:",
        paste(colnames(Y_tmp), collapse = ", "),
        ", instead.\n"
      ))
    }

    # embed: get lagged values
    X_train <- embed(Y_tmp, dimension = p + 1)[, -(1:M)]
    lag_idx <- rep(1:p, each = M)
    var_idx <- rep(colnames(Y_tmp), times = p)
    colnames(X_train) <- paste0(var_idx, ".t-", lag_idx)
    
    Y_train <- Y_tmp[-c(1:p), ]
    TT <- Traw - p
    Y_fit_BART <- Y_train * 0
    
    colnames(Y_train) <- colnames(Y_tmp)
    
    datamat <- data.frame(cbind(Y_train, X_train))
    has_test <- FALSE
    
    variable_names <- colnames(Y_train)
    index_names <- rownames(Y_train)
    
    if (TT != nrow(Y_train) | TT != nrow(X_train)) {
      stop("Something went wrong: Tobs != nrow(Y). \n")
    }

    ###--------------------------------------------------------------------------###
    ###-------------------------- MLE estimates ---------------------------------###
    ###--------------------------------------------------------------------------###
    
    XtX <- t(X_train) %*% X_train
    XtY <- t(X_train) %*% Y_train
    
    # Initialization of the covariance dependes on the high dimensionality
    if (TT <= K) {
      beta_ols <- MASS::ginv(XtX) %*% XtY
      y_hat_mle <- X_train %*% beta_ols
      resid_mle <- Y_train - y_hat_mle
      sigma2hat <- (t(resid_mle) %*% (resid_mle)) / TT
    } else {
      beta_ols <- solve(XtX) %*% XtY
      y_hat_mle <- X_train %*% beta_ols
      resid_mle <- Y_train - X_train %*% beta_ols
      sigma2hat <- (t(resid_mle) %*% (resid_mle)) / TT
    }

    if(variance_prior == "csv"){
        U <- diag(M)
    }

    ###--------------------------------------------------------------------------###
    ###-------------------- Preprocessing the BART objects  ---------------------###
    ###--------------------------------------------------------------------------###
    
    # --- General params ---
    
    general_params_default <- list(
      cutpoint_grid_size = 100,
      standardize = TRUE,
      quantile_transform = FALSE,
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
      num_threads = -1,
      k_hyper = 2.,
      a_dart = 0.5,
      b_dart = 1.0,
      lambda_1 = 1.,
      lambda_2 = 0.5
    )
    general_params_updated <- preprocessParams(general_params_default, general_params)
    
    # --- Mean-forest  ---
    mean_forest_params_default <- list(
      num_trees = 200,
      alpha = 0.95,
      beta = 2.0,
      min_samples_leaf = 5,
      max_depth = 10,
      fixed_sigma2_leaf = FALSE,
      sample_sigma2_leaf = TRUE,     # τ sampling 
      sigma2_leaf_init = NULL,        # init τ
      sigma2_leaf_shape = 3,       
      sigma2_leaf_scale = NULL,    
      keep_vars = NULL,
      drop_vars = NULL,
      num_features_subsample = NULL
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
    quantile_transform   <- general_params_updated$quantile_transform
    k_hyper              <- general_params_updated$k_hyper
    
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
    fixed_sigma2_leaf    <- mean_forest_params_updated$fixed_sigma2_leaf
    

    ###-------------------- Data preprocessing --------------------###
    train_cov_preprocess_list <- preprocessTrainData(X_train)
    X_train_metadata <- train_cov_preprocess_list$metadata
    X_train <- train_cov_preprocess_list$data
    original_var_indices <- X_train_metadata$original_var_indices
    feature_types <- X_train_metadata$feature_types

    var_y <-diag(var(Y_train))

    if (is.null(sigma2_init)) {
      sigma2_init <- 1.0 * var_y
    }
    
    if (is.null(b_leaf)) {
      b_leaf <- var_y / (2 * num_trees_mean)
    }
    
    if (is.null(sigma2_leaf_init)) {
      sigma2_leaf_init <- as.matrix(
        2 * var_y / (num_trees_mean)
      )
      current_leaf_scale <- sigma2_leaf_init
    } else {
      current_leaf_scale <- sigma2_leaf_init
    }
    
    current_sigma2 <- sigma2_init
    
    # Random number generator (std::mt19937)
    if (is.null(random_seed)) {
      random_seed = sample(1:10000, 1, FALSE)
    }
    rng <- stochtree::createCppRNG(random_seed)
    feature_types <- as.integer(feature_types)

    if (fixed_sigma2_leaf) {
      sample_sigma2_leaf <- FALSE

      ymin <- apply(Y_train, 2, min)
      ymax <- apply(Y_train, 2, max)

      node_scale <- (ymax - ymin) / 2       
      sigma_mu   <- node_scale / (k_hyper * sqrt(num_trees_mean))
      current_leaf_scale <- sigma_mu^2      
    }

    if(mean_prior == "sparse"){
        a_dart <- mean_forest_params_updated$a_dart
        b_dart <- mean_forest_params_updated$b_dart
        
        alpha_dart <-1.
        rho_dart <- ncol(X_train)
    }

    if (mean_prior == "minnesota"){
        lambda_1 <- mean_forest_params_updated$lambda_1
        lambda_2 <- mean_forest_params_updated$lambda_2
        lag_index <- matrix(0,M,p) # For the minessota prior.  
        
        for(i in 1:M){
            lag_index[i,] <- seq(i, K, by = M)
        }
    }


    if (quantile_transform) {
      ecdfs <- apply(X_train, 2, ecdf)  # list of ecdf functions
      for (j in seq_along(ecdfs)) {
        f <- ecdfs[[j]]
        X_train[, j] <- f(X_train[, j])
      }
    }
    
    # Unpack model type information
    leaf_model_mean_forest = 0
    leaf_dimension = 1
    is_leaf_constant = TRUE
    leaf_regression = FALSE
    init_values_mean_forest <- 0.
        
    ###-------------------------------Stochtree Objects--------------------------###
    
    forest_dataset_train_ls <- c()
    forest_model_config_mean_ls <- c()
    active_forest_mean_ls <- c()
    forest_samples_mean_ls <- c()
    global_model_config_ls <- c()
    forest_model_mean_ls <- c()
    variable_weights_mean_ls <- matrix(0, nrow = M, ncol = K)
    
    for (mm in 1:M) {
        variable_weights_mean_ls[mm, ] <- rep(1 / K, K)
        forest_dataset_train_ls <- c(
        forest_dataset_train_ls,
        createForestDataset(X_train, basis = NULL, variance_weights = rep(1, TT))
        )
        forest_model_config_mean_ls <- c(
        forest_model_config_mean_ls,
        stochtree::createForestModelConfig(
            feature_types = feature_types,
            num_trees = num_trees_mean,
            num_features = ncol(X_train),
            num_observations = nrow(X_train),
            variable_weights = variable_weights_mean_ls[mm, ], # this actually not a list is a matrix.
            leaf_dimension = leaf_dimension,
            alpha = alpha_mean,
            beta = beta_mean,
            min_samples_leaf = min_samples_leaf_mean,
            max_depth = max_depth_mean,
            leaf_model_type = leaf_model_mean_forest,
            leaf_model_scale = current_leaf_scale[mm],
            cutpoint_grid_size = cutpoint_grid_size,
            num_features_subsample = NULL
        )
        )
        forest_samples_mean_ls <- c(
        forest_samples_mean_ls,
        stochtree::createForestSamples(
            num_trees_mean,
            leaf_dimension,
            is_leaf_constant,
            FALSE
        )
        )
        active_forest_mean_ls <- c(
        active_forest_mean_ls,
        stochtree::createForest(
            num_trees_mean,
            leaf_dimension,
            is_leaf_constant,
            FALSE
        )
        )
        
        if(variance_prior != "const" ){
        global_model_config_ls <- c(
            global_model_config_ls,
            stochtree::createGlobalModelConfig(1.))
        }else{
        global_model_config_ls <- c(
            global_model_config_ls,
            stochtree::createGlobalModelConfig(current_sigma2[mm]))
        }
        
        forest_model_mean_ls <- c(
        forest_model_mean_ls,
        stochtree::createForestModel(
            forest_dataset_train_ls[[mm]],
            forest_model_config_mean_ls[[mm]],
            global_model_config_ls[[mm]]
        )
        )
    }
    
    # Initialize the leaves of each tree in the mean forest, for each equation
    for (mm in 1:M) {
      outcome_train <- stochtree::createOutcome(Y_train[, mm])
      forest_dataset_train <- forest_dataset_train_ls[[mm]]
      active_forest_mean <- active_forest_mean_ls[[mm]]
      forest_model_mean <- forest_model_mean_ls[[mm]]

      active_forest_mean$prepare_for_sampler(
        forest_dataset_train,
        outcome_train,
        forest_model_mean,
        leaf_model_mean_forest,
        init_values_mean_forest
      )

      active_forest_mean_ls[[mm]] <- active_forest_mean
    }

    ###--------------------------------------------------------------------------###
    ###-------------------------- Covariance OBJECTS ----------------------------###
    ###--------------------------------------------------------------------------###


    if (variance_prior == "fsv" || variance_prior == "csv") {
        sv_priors <- list()
        sv_params <- list()
        h_latent <- list()

        for (mm in 1:M) {
          sv_params[[mm]] <- list(
            # mu      = -3 + rnorm(1),
            mu      = 0 + rnorm(1, sd = 0.1),
            phi     = 0.8 + pmin(rnorm(1, sd = 0.06), 0.095),
            sigma   = 0.1 + rgamma(1, shape = 1, rate = 10),
            nu      = Inf,
            rho     = 0,
            beta    = NA,
            latent0 = 0
          )
          sv_params[[mm]]$latent0 <- sv_params[[mm]]$mu + rnorm(1)
          h_latent[[mm]] <- sv_params[[mm]]$mu + rnorm(TT)
          sv_priors[[mm]] <- stochvol::specify_priors(
            mu     = stochvol::sv_normal(mean = 0, sd = 1),
            phi    = stochvol::sv_beta(shape1 = 22, shape2 = 1.5),
            sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 50),
            # sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 1/(2*0.01)),
            nu     = stochvol::sv_infinity(),
            rho    = stochvol::sv_constant(0)
          )
        }
        current_sigma2 <- rep(1., M)
    }

    if(variance_prior == "fsv"){
        # Number of factors
        if(is.null(num_factors)){
        Q <- upper_bound(M)
        
        if (Q == 0 || Q == 1) {
          Q <- 1
        } # at least 1 factors
        }else{
          Q <- num_factors
        }
        Lambda <- matrix(0, M, Q)
        Ft     <- matrix(0, TT, Q)
        # Lambda <- matrix(rnorm(M * Q, sd = 0.3), nrow = M, ncol = Q)
        # Ft     <- matrix(rnorm(TT * Q, sd = 0.1), nrow = TT, ncol = Q)
        
        Omega <- matrix(0, TT, Q)
        theta_Lambda <- matrix(1,M, Q)
        
        fsv_priors <- list()
        fsv_params <- list()
        fsv_h_latent <- list()
        
        for (qq in 1:Q) {
          fsv_params[[qq]] <- list(
            mu      = 0,
            phi     = 0.8 + pmin(rnorm(1, sd = 0.06), 0.095),
            sigma   = 0.1 + rgamma(1, shape = 1, rate = 10),
            nu      = Inf,
            rho     = 0,
            beta    = NA,
            latent0 = rnorm(1)
          )
          fsv_h_latent[[qq]] <- fsv_params[[qq]]$mu + rnorm(TT)
          fsv_priors[[qq]] <- stochvol::specify_priors(
            mu     = stochvol::sv_constant(0),
            phi    = stochvol::sv_beta(shape1 = 22, shape2 = 1.5),
            # sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 0.5),
            sigma2 = stochvol::sv_gamma(shape = 0.5, rate = 50),
            nu     = stochvol::sv_infinity(),
            rho    = stochvol::sv_constant(0)
          )
        }
        
        current_sigma2 <- rep(1., M)
        # HS prior for factor loadings (by column)
        if(cov_shrinkage == "hs"){
            lambda_Lambda <- matrix(1, M, Q)
            nu_Lambda <- matrix(1, M, Q)
            zeta_Lambda <- rep(1, Q)
            tau_Lambda <- rep(1, Q)
        }else if(cov_shrinkage == "ng"){
            theta_Lambda <- matrix(1, M, Q)   # prior variances for Λ
            lambda_ng    <- rep(1, Q)         # global column params for NG
            aj <- rep(0.1,   Q)
            cj <- rep(0.001, Q)
            dj <- rep(0.001, Q)
        }

    }

    if (variance_prior == "csv") {
        # Lower-triangular contemporaneous matrix A with diag(A)=1 such that:
        #   A * eta_t = eps_t,   eps_t ~ N(0, diag(exp(h_t)))
        # and eta_t = y_t - f_t.
        A <- diag(M)

        # Prior-variance containers for A off-diagonals (row-wise)
        theta_A <- matrix(1, M, M)

        if (cov_shrinkage == "hs") {
            lambda_A <- matrix(1, M, M)  # local scales
            nu_A     <- matrix(1, M, M)  # aux variables
            tau_A    <- rep(1, M)        # row-global scales
            zeta_A   <- rep(1, M)        # aux variables
        } else if (cov_shrinkage == "ng") {
            lambda_ng_A <- rep(1, M)     # row-global NG parameter (optional book-keeping)
            aj_A <- rep(0.1, M)
            cj_A <- rep(0.001, M)
            dj_A <- rep(0.001, M)
        }
    }
   
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
    
    ###--------------------------------------------------------------------------###
    ###-------------------- Storage ---------------------------------------------###
    ###--------------------------------------------------------------------------###
    
    y_hat_store <- array(NA, dim = c(num_saved, TT, M))
    sigma_store <- array(NA, dim = c(num_saved, TT, M))
    var_count_matrix <- array(NA, dim = c(num_saved, K, M))
    
    if (variance_prior == "fsv") {
      Lambda_store <- array(NA, dim = c(num_saved, M, Q))
      Ot_store <- array(NA, dim = c(num_saved, TT, Q))
      Factor_store <- array(NA, dim = c(num_saved, TT, Q))
      sv_params_mcmc <- array(NA, dim = c(num_saved, 3, M))
      sv_params_mat <- matrix(NA, M, 4)
      colnames(sv_params_mat) = c("mu", "phi", "sigma", "h")

      fsv_params_mcmc <- array(NA, dim = c(num_saved, 3, Q))
      fsv_params_mat <- matrix(NA, Q, 4)
      colnames(fsv_params_mat) <- c("mu", "phi", "sigma", "h")

      Ht_store <- array(NA, dim = c(num_saved, TT, M))
      H_t <- matrix(0, TT, M)
    }
    
    if (variance_prior == "csv") {
        A_store <- array(NA, dim = c(num_saved, M, M))
        sv_params_mcmc <- array(NA, dim = c(num_saved, 3, M))
        sv_params_mat <- matrix(NA, M, 4)
        colnames(sv_params_mat) = c("mu", "phi", "sigma", "h")

        Ht_store <- array(NA, dim = c(num_saved, TT, M))
        H_t <- matrix(0, TT, M)
    }
    
    global_var_samples <- matrix(0, num_saved, M)
    colnames(global_var_samples) <- variable_names
    
    if (!is.null(Y_test)) {
      if (is.null(n_ahead)) {
        n_ahead <- 12
      }
      predictions <- array(
        as.numeric(NA),
        c(num_saved, n_ahead, M),
        dimnames = list(
          paste0("mcmc_", 1:num_saved),
          paste0("t+", 1:n_ahead),
          variable_names
        )
      )
      sigma_predictions <- array(NA, c(num_saved, n_ahead, M))
      LPL_draws <- matrix(as.numeric(NA), num_saved, n_ahead)
      colnames(LPL_draws) <- paste0("t+", 1:n_ahead)
      PL_univariate_draws <- array(
        as.numeric(NA),
        c(num_saved, n_ahead, M),
        dimnames = list(
          paste0("mcmc_", 1:num_saved),
          paste0("t+", 1:n_ahead),
          variable_names
        )
      )

      has_test <- TRUE
    }


    if(num_warmstart > 0){
      cat("Running warmstart iterations... \n")
      for(ws in 1:num_warmstart){
        for(mm in 1:M){
          outcome_train <- stochtree::createOutcome(Y_train[, mm] - Y_fit_BART[, mm])
          forest_dataset_train <- forest_dataset_train_ls[[mm]]
          active_forest_mean <- active_forest_mean_ls[[mm]]
          forest_model_mean <- forest_model_mean_ls[[mm]]
          forest_samples_mean <- forest_samples_mean_ls[[mm]]
          forest_model_config_mean <- forest_model_config_mean_ls[[mm]]
          global_model_config <- global_model_config_ls[[mm]]
          
          forest_model_mean$sample_one_iteration(
            forest_dataset = forest_dataset_train,
            residual = outcome_train,
            forest_samples = forest_samples_mean,
            active_forest = active_forest_mean,
            rng = rng,
            forest_model_config = forest_model_config_mean,
            global_model_config = global_model_config,
            num_threads = num_threads,
            keep_forest = FALSE,
            gfr = TRUE
          )
          current_sigma2[mm] <- stochtree::sampleGlobalErrorVarianceOneIteration(
                outcome_train,
                forest_dataset_train,
                rng,
                a_global,
                b_global
                )
                global_model_config_ls[[
                mm
                ]]$update_global_error_variance(current_sigma2[mm])
        
        leaf_scale_double <- sampleLeafVarianceOneIteration(active_forest_mean,rng,
                                                                    a_leaf,b_leaf[mm])
        current_leaf_scale[mm] <- leaf_scale_double
        forest_model_config_mean_ls[[mm]]$update_leaf_model_scale(as.matrix(current_leaf_scale[mm]))
        #Updating
        forest_model_mean_ls[[mm]] <- forest_model_mean  
        Y_fit_BART[, mm] <- active_forest_mean$predict(forest_dataset_train)
        }
      }
      cat("Warmstart completed. \n")
    }
    
    ###--------------------------------------------------------------------------###
    ###----------------------------- MCMC sampler -------------------------------###
    ###----------------------------- for VARBART  -------------------------------###
    ###--------------------------------------------------------------------------###

    pb = txtProgressBar(min = 0, max = num_samples, style = 3)
    start = Sys.time()
    for (i in 1:num_samples) {    
        index_saved <- if (i %in% save_set) match(i, save_set) else NA_integer_
        
        is_mcmc <- i > num_burnin
        
        if (is_mcmc) {
            mcmc_counter <- i - num_burnin
            if (has_test) {
                keep_sample <- FALSE
            } else {
                keep_sample <- (mcmc_counter %% num_thin == 0)
            }
            } else {
            if (has_test || !keep_burnin) {
                keep_sample <- FALSE
            } else {
                keep_sample <- TRUE
            }
        }
        
        for (mm in 1:M) {
            # ------------------ Construct equation-specific response ------------------ #
            if (variance_prior == "fsv") {
                # remove factor component from equation mm
                y_star <- Y_train[, mm] - drop(tcrossprod(Ft, Lambda[mm, , drop = FALSE]))
            } else if (variance_prior == "csv") {
                # Carriero-style triangularization: y*_m = y_m + sum_{j<m} a_{m,j} * eta_j
                if (mm == 1) {
                    y_star <- Y_train[, mm]
                } else {
                    eta_prev <- Y_train[, 1:(mm-1), drop = FALSE] - Y_fit_BART[, 1:(mm-1), drop = FALSE]
                    y_star <- as.numeric(Y_train[, mm] + eta_prev %*% matrix(A[mm, 1:(mm-1)], ncol = 1))
                }
            } else {
                y_star <- Y_train[, mm]
            }

            outcome_train <- stochtree::createOutcome(y_star - Y_fit_BART[, mm])

            # # Sampling Model Mean Function:
            # if (variance_prior == "fsv") {
            #     Y_ <- Y_train - tcrossprod(Ft, Lambda)
            # } else {
            #     Y_ <- Y_train
            # }

            # outcome_train <- stochtree::createOutcome(Y_[, mm] - Y_fit_BART[, mm])

            ###------------------------------ Step 1: -----------------------------------###
            ###---------------------- Sample BART objects--------------------------------###
            ###--------------------------------------------------------------------------###

            #Unpacking the objects for each equation:

            forest_dataset_train <- forest_dataset_train_ls[[mm]]
            active_forest_mean <- active_forest_mean_ls[[mm]]
            forest_model_mean <- forest_model_mean_ls[[mm]]
            forest_samples_mean <- forest_samples_mean_ls[[mm]]
            forest_model_config_mean <- forest_model_config_mean_ls[[mm]]
            global_model_config <- global_model_config_ls[[mm]]

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
            variable_count_splits <- active_forest_mean$get_forest_split_counts(ncol(
                X_train
            ))
            Y_fit_BART[, mm] <- active_forest_mean$predict(forest_dataset_train)

        
            # Global Error Variance
            if (variance_prior == "const") {
                current_sigma2[mm] <- stochtree::sampleGlobalErrorVarianceOneIteration(
                outcome_train,
                forest_dataset_train,
                rng,
                a_global,
                b_global
                )
                global_model_config_ls[[
                mm
                ]]$update_global_error_variance(current_sigma2[mm])
            }

            # Leaf Node Variance
            if (sample_sigma2_leaf) {
                leaf_scale_double <- sampleLeafVarianceOneIteration(active_forest_mean,rng,
                                                                    a_leaf,b_leaf[mm])
                current_leaf_scale[mm] <- leaf_scale_double
                forest_model_config_mean_ls[[mm]]$update_leaf_model_scale(as.matrix(current_leaf_scale[mm]))
            }

            if (i > floor(num_burnin / 2)) {
                if(mean_prior == "sparse"){
            
                    draw_dart <- sample_dart_splits_one_iteration(variable_count_splits, alpha_dart, rng)
                    log_prob_vector <- draw_dart$lpv
                    variable_prob_splits <- exp(log_prob_vector)

                    alpha_sampler <- sample_alpha_one_iteration(log_prob_vector, a_dart, b_dart, rho_dart, rng)
                    alpha_dart    <- alpha_sampler$alpha
            
                }
                if(mean_prior == "minnesota"){
                    log_probability_spits <- draw_minessota_split(variable_count_splits, mm,lag_index ,diag(sigma2hat),lambda_1, lambda_2, rng)
                    variable_prob_splits      <- exp(log_probability_spits)
                }
                #Updating variable weights
                forest_model_config_mean_ls[[mm]]$update_variable_weights(variable_prob_splits)
            }
            #Updating
            forest_model_mean_ls[[mm]] <- forest_model_mean

            if (!is.na(index_saved)) {
                global_var_samples[index_saved, mm] <- current_sigma2[mm] * (Ysd[mm]^2)
                var_count_matrix[index_saved, , mm] <- variable_count_splits
            }
        } # end 1:M BART loop
        
        if (variance_prior == "csv") {
          # reduced-form residuals
          eta <- Y_train - Y_fit_BART  # TT x M

          # structural shocks eps = A * eta  (A is lower-triangular with diag=1)
          eps <- eta
          for (m in 2:M) {
            eps[, m] <- as.numeric(eta[, m] + eta[, 1:(m-1), drop = FALSE] %*% matrix(A[m, 1:(m-1)], ncol = 1))
          }

          # ---- update each row m=2..M of A via weighted regression: eta_m = -eta_<m a_m + eps_m ----
          for (m in 2:M) {
            w <- as.numeric(exp(-0.5 * h_latent[[m]]))  # TT
            X <- (-eta[, 1:(m-1), drop = FALSE]) * w    # TT x (m-1)
            y <- eta[, m] * w                            # TT

            # prior variances for this row's coefficients
            if (cov_shrinkage == "hs") {
              theta_A[m, 1:(m-1)] <- (tau_A[m]^2) * (lambda_A[m, 1:(m-1)]^2)
            }
            V0 <- as.numeric(theta_A[m, 1:(m-1), drop = TRUE])

            S <- crossprod(X)
            diag(S) <- diag(S) + 1 / V0
            R <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-8, ncol(X))))

            b <- crossprod(X, y)
            mu <- backsolve(R, forwardsolve(t(R), b))
            draw <- mu + backsolve(R, rnorm(ncol(X)))

            A[m, 1:(m-1)] <- as.numeric(draw)

            # ---- shrinkage hyperparameter updates for A ----
            if (cov_shrinkage == "hs") {
              hs <- sample_hs(
                bdraw = as.numeric(A[m, 1:(m-1)]),
                lambda.hs = as.numeric(lambda_A[m, 1:(m-1)]),
                nu.hs = as.numeric(nu_A[m, 1:(m-1)]),
                tau.hs = tau_A[m],
                zeta.hs = zeta_A[m]
              )
              lambda_A[m, 1:(m-1)] <- hs$lambda
              nu_A[m, 1:(m-1)]     <- hs$nu
              tau_A[m]             <- hs$tau
              zeta_A[m]            <- hs$zeta
              theta_A[m, 1:(m-1)]  <- (tau_A[m]^2) * (lambda_A[m, 1:(m-1)]^2)
            } else if (cov_shrinkage == "ng") {
              upd <- sample_ng_col(
                facloads_j = as.numeric(A[m, 1:(m-1)]),
                tau2_j     = as.numeric(theta_A[m, 1:(m-1)]),
                aj = aj_A[m],
                cj = cj_A[m],
                dj = dj_A[m]
              )
              lambda_ng_A[m]       <- upd$lambda_j
              theta_A[m, 1:(m-1)]  <- upd$tau2_j
            }
          }

          theta_A[theta_A < 1e-8] <- 1e-8
        }

        
        ###--------------------- Step 2: Idiossincratic Vol -------------------------###
        ###--------------- Sample stochastic volatility parameters ------------------###
        ###--------------------------------------------------------------------------###
        
        if (variance_prior == "fsv" || variance_prior == "csv") {
            
            if(variance_prior == "fsv"){
                shocks = Y_train - tcrossprod(Ft, Lambda) - Y_fit_BART
            }else if(variance_prior == "csv"){
                eta <- Y_train - Y_fit_BART
                shocks <- eta
                for (m in 2:M) {
                    shocks[, m] <- as.numeric(eta[, m] + eta[, 1:(m-1), drop = FALSE] %*% matrix(A[m, 1:(m-1)], ncol = 1))
                }
            }else{
                shocks = Y_train - Y_fit_BART
            }
            for (mm in 1:M) {
                svdraw_mm = stochvol::svsample_general_cpp(
                shocks[, mm],
                startpara = sv_params[[mm]],
                startlatent = h_latent[[mm]],
                priorspec = sv_priors[[mm]]
                )
                sv_params[[mm]][c("mu", "phi", "sigma")] = as.list(svdraw_mm$para[, c(
                "mu",
                "phi",
                "sigma"
                )])

                sv_params_mat[mm, ] = c(
                svdraw_mm$para[, c("mu", "phi", "sigma")],
                svdraw_mm$latent[TT]
                )

                h_raw <- as.numeric(svdraw_mm$latent)              # log var that stochvol fits
                mu_h  <- as.numeric(svdraw_mm$para[, "mu"])        # stationary mean of h_raw

                current_sigma2[mm] <- exp(mu_h)                     # σ^2 = exp(μ)
                u_t <- h_raw - mu_h                             # mean-zero log-vol
                u_t <- pmax(u_t, log(1e-6))                     
                wts <- exp(u_t)                                 # variance weights
                h_t <- h_raw

                if (i %in% save_set) {
                    sv_params_mcmc[index_saved, , mm] = svdraw_mm$para[, c(
                        "mu",
                        "phi",
                        "sigma"
                    )]
                }

                forest_dataset_train_ls[[mm]]$update_variance_weights(variance_weights = wts)
                global_model_config_ls[[mm]]$update_global_error_variance(current_sigma2[mm])

                h_latent[[mm]] <- h_raw
                H_t[, mm] <- h_raw
            }
            H_t[H_t < log(1e-6)] <- log(1e-6)
        }
        
        ###------------------------------ Step 3: -----------------------------------###
        ###---------------- Sample the factor loadings and the factors --------------###
        ###--------------------------------------------------------------------------###
        
        ###--------- Step 3.4: Sample Factors Idiossincratic Volatilities------------###
        if (variance_prior == "fsv") {
            ## ------------------- Step 3.1: Sample Factors ------------------- ##
            resid_factors <- Y_train - Y_fit_BART
            Ft <- sample_factors(resid_factors, Lambda, H_t, Omega) # T x Q

            ## ------------------- Step 3.2: Sample Loadings ------------------- ##

            Lambda <- sample_loadings(resid_factors, Ft, H_t, theta_Lambda)

            ## --------- Step 3.3: Columnwise shrinkage on Λ --------- ##
            for (qq in 1:Q) {
              if(cov_shrinkage == "hs"){
                hs_draw <- sample_hs(
                bdraw = as.numeric(Lambda[, qq]),
                lambda.hs = as.numeric(lambda_Lambda[, qq]),
                nu.hs = nu_Lambda[, qq],
                tau.hs = tau_Lambda[qq],
                zeta.hs = zeta_Lambda[qq]
                )
                theta_Lambda[, qq] <- hs_draw$psi
                lambda_Lambda[, qq] <- hs_draw$lambda
                nu_Lambda[, qq] <- hs_draw$nu
                tau_Lambda[qq] <- hs_draw$tau
                zeta_Lambda[qq] <- hs_draw$zeta
              }else if(cov_shrinkage == "ng"){
                upd <- sample_ng_col(facloads_j = as.numeric(Lambda[, qq]),
                                      tau2_j     = theta_Lambda[, qq], 
                                      aj = aj[qq],
                                      cj = cj[qq],
                                      dj = dj[qq])
                lambda_ng[qq]      <- upd$lambda_j   # column-global λ_j
                theta_Lambda[, qq] <- upd$tau2_j     # updated τ_{ij}^2
              }
            }
            theta_Lambda[theta_Lambda < 1e-5] <- 1e-5

            ## --------------- Step 3.4: Sample factor SV states Ω -------------- ##

            for (qq in 1:Q) {
                fsvdraw_qq <- stochvol::svsample_general_cpp(
                Ft[, qq],
                startpara = fsv_params[[qq]],
                startlatent = fsv_h_latent[[qq]],
                priorspec = fsv_priors[[qq]]
                )
                fsv_params[[qq]][c("mu", "phi", "sigma")] <- as.list(fsvdraw_qq$para[, c(
                "mu",
                "phi",
                "sigma"
                )])
                fsv_h_latent[[qq]] <- fsvdraw_qq$latent
                fsv_params_mat[qq, ] <- c(
                fsvdraw_qq$para[, c("mu", "phi", "sigma")],
                fsvdraw_qq$latent[TT]
                )
                if (i %in% save_set) {
                fsv_params_mcmc[index_saved, , qq] <- fsvdraw_qq$para[, c(
                    "mu",
                    "phi",
                    "sigma"
                )]
                }
                Omega[, qq] <- fsvdraw_qq$latent # (log-variance path)
            }
            Omega[Omega < log(1e-6)] <- log(1e-6)
        } # end FSV
        

        if (!is.na(index_saved)){

            y_hat_store[index_saved, , ] <- (Y_fit_BART * t(matrix(Ysd, M, TT))) +
                t(matrix(Ymu, M, TT))
            if (variance_prior == "const") {
                sigma_store[index_saved, , ] <- current_sigma2 * (Ysd^2)
            }
            if (variance_prior == "fsv") {
                # Ht_store[index_saved, , ] <- exp(H_t * 0.5) * t(matrix(Ysd, M, TT))
                Ht_store[index_saved, , ]     <- exp(H_t) 
                Ot_store[index_saved, , ]     <- exp(Omega)
                Lambda_store[index_saved, , ] <- Lambda
                Factor_store[index_saved,,]    <- Ft
            }
            if(variance_prior == "csv"){
                Ht_store[index_saved, , ] <- exp(H_t)
                A_store[index_saved, , ] <- A 
            }
        
        if (has_test) {
            ###------------------------------------------------------------------------------###
            ###----------------------------- Forecasting ------------------------------------###
            ###------------------------------------------------------------------------------###
            X_hat <- datamat[nrow(datamat), 1:K]
            if (variance_prior == "fsv") {
                HT <- H_t[TT, ]
                OT <- Omega[TT, ]
            } else {
                HT <- current_sigma2
            }

            if (variance_prior == "csv") {
                HT <- H_t[TT, ]
            }
            
            X_fore_k <- as.matrix(X_hat)
            forest_dataset <- stochtree::createForestDataset(X_fore_k)
            mean_forecast <- matrix(0, M)
            
            for (k in 1:n_ahead) {
            # 1-step-ahead mean from BART
            for (mm in 1:M) {
                active_forest_sample_mean <- active_forest_mean_ls[[mm]]
                mean_forecast[mm] <- active_forest_sample_mean$predict_raw(forest_dataset)
            }
            
            if (variance_prior == "fsv") {
                # AR(1) updates for log-vols: HT (idiosyncratic), OT (factor)
                HT <- (sv_params_mat[, 1] +
                        sv_params_mat[, 2] * (HT - sv_params_mat[, 1]) +
                        sv_params_mat[, 3] * rnorm(M))
                
                OT <- (fsv_params_mat[, 1] +
                        fsv_params_mat[, 2] * (OT - fsv_params_mat[, 1]) +
                        fsv_params_mat[, 3] * rnorm(Q))
                
                # Sigma in standardized space
                Sigma_fore <- Lambda %*% diag(exp(OT)) %*% t(Lambda) + diag(exp(HT))
            }else if(variance_prior == "csv"){
                HT <- (sv_params_mat[, 1] +
                        sv_params_mat[, 2] * (HT - sv_params_mat[, 1]) +
                        sv_params_mat[, 3] * rnorm(M))
                # Sigma in standardized space
                Ainv <- forwardsolve(A, diag(M))
                Sigma_fore <- Ainv %*% diag(exp(HT)) %*% t(Ainv)
            }else {
                h_fore <- HT # should be variances
                Sigma_fore <- diag(h_fore)
            }

            

            
            ## ---------- Draw forecast in original scale (unchanged) ----------
            predictions[index_saved, k, ] <- MASS::mvrnorm(
                1,
                mean_forecast,
                Sigma_fore
            )
            
            sigma_predictions[index_saved, k, ] <- diag(Sigma_fore)
            
            Y_obs_vec <- as.vector(Y_test[k, ])
            mean_forecast_scaled <- (mean_forecast * Ysd) + Ymu
            
            if (variance_prior == "fsv") {
                ## ---------- Woodbury-based predictive log-likelihood ----------
                # Work in standardized space
                y_std  <- (Y_obs_vec - Ymu) / Ysd
                mu_std <- as.vector(mean_forecast) # mean_forecast is already standardized in your code
                
                # Calculate density on STANDARDIZED scale
                lpl_std <- logdmvnorm_fsv_woodbury(
                y      = y_std,
                mu     = mu_std,
                Lambda = Lambda,
                h_idio = HT,
                h_fac  = OT
                )
                
                
                LPL_draws[index_saved, k] <- lpl_std - sum(log(Ysd))
                
                ## Univariate draws correction
                # Calculate variance in standardized space
                D    <- exp(HT)
                fvar <- exp(OT)
                var_std <- D + rowSums((Lambda^2) %*% diag(fvar))
                
                # Convert to Original Scale Variance
                var_orig <- (Ysd^2) * var_std
                
                PL_univariate_draws[index_saved, k, ] <- stats::dnorm(
                Y_obs_vec,
                mean_forecast_scaled,
                sqrt(var_orig),
                log = TRUE
                )
                
            }else if(variance_prior == "csv"){
                
                ## ---------- Predictive log-likelihood ----------
                # Work in standardized space
                y_std  <- (Y_obs_vec - Ymu) / Ysd
                mu_std <- as.vector(mean_forecast) # mean_forecast is already standardized in your code
                Sigma_std <- Ainv %*% diag(exp(HT)) %*% t(Ainv)
                # Calculate density on STANDARDIZED scale
                lpl_std <- mvtnorm::dmvnorm(
                y_std,
                mu_std,
                Sigma_std,
                log = TRUE
                )
                
                LPL_draws[index_saved, k] <- lpl_std - sum(log(Ysd))
                
                ## Univariate draws correction
                # Calculate variance in standardized space
                var_std <- diag(Sigma_std)
                
                # Convert to Original Scale Variance
                var_orig <- (Ysd^2) * var_std
                
                PL_univariate_draws[index_saved, k, ] <- stats::dnorm(
                Y_obs_vec,
                mean_forecast_scaled,
                sqrt(var_orig),
                log = TRUE
                )
            }else{   
                Sigma_fore_scaled <- diag(Ysd) %*% Sigma_fore %*% diag(Ysd)
                
                LPL_draws[index_saved, k] <- mvtnorm::dmvnorm(
                Y_obs_vec,
                mean_forecast_scaled,
                Sigma_fore_scaled,
                log = TRUE
                )
                
                PL_univariate_draws[index_saved, k, ] <- stats::dnorm(
                Y_obs_vec,
                mean_forecast_scaled,
                sqrt(diag(Sigma_fore_scaled)),
                log = TRUE
                )
            }
            
            # ---------- Update regressors for multi-step forecasts ----------
            if (k < n_ahead) {
                if (p == 1) {
                X_fore_k[1, ] <- predictions[index_saved, k, ]
                } else {
                X_fore_k[1, ] <- c(
                    predictions[index_saved, k, ],
                    X_fore_k[1:((p - 1) * M)]
                )
                }
                forest_dataset <- stochtree::createForestDataset(X_fore_k)
                }
                } # end 1:k
            } # has test
        } # end if index_saved
        
        iter_update <- 250
        setTxtProgressBar(pb, i)
        if (i %% iter_update == 0) {
        end <- Sys.time()
        message(paste0(
            "\n Average time for single draw over last ",
            iter_update,
            " draws ",
            round(as.numeric(end - start) / iter_update, digits = 4),
            " seconds, currently at draw ",
            i
        ))
        start <- Sys.time()
        }
    } # end 1:num_samples
    close(pb)
    
    dimnames(y_hat_store) <- list(
        paste0("mcmc_", 1:num_saved),
        index_names,
        variable_names
    )
    dimnames(var_count_matrix) <- list(
        paste0("mcmc_", 1:num_saved),
        colnames(X_train),
        variable_names
    )
    
    model_params <- list(
        "sigma2_init" = sigma2_init,
        "a_global" = a_global,
        "b_global" = b_global,
        "a_leaf" = a_leaf,
        "b_leaf" = b_leaf,
        "outcome_mean" = Ymu,
        "outcome_scale" = Ysd,
        "leaf_dimension" = leaf_dimension,
        "is_leaf_constant" = is_leaf_constant,
        "leaf_regression" = leaf_regression,
        "num_covariates" = ncol(X_train),
        "num_burnin" = num_burnin,
        "num_mcmc" = num_mcmc,
        "sample_sigma2_global" = sample_sigma2_global,
        "variance_prior" = variance_prior
    )
    
    result <- list(
        "model_params" = model_params,
        "train_set_metadata" = X_train_metadata,
        "var_count" = var_count_matrix,
        "lags" = p,
        "M" = M,
        "K" = K,
        "Y_raw" = data
    )
    
    result[["y_hat_train"]] = y_hat_store
    
    if (variance_prior == "fsv") {
        result[["Lambda_store"]] = Lambda_store
        result[["Factor_store"]] = Factor_store
        result[["Ot_store"]] = Ot_store
        result[["Ht_store"]] = Ht_store
        result[["sv_params_mcmc"]] = sv_params_mcmc
        result[["fsv_params_mcmc"]] = fsv_params_mcmc
    }
    if (variance_prior == "csv") {
        result[["Ht_store"]] = Ht_store
        result[["sv_params_mcmc"]] = sv_params_mcmc
        result[["A_store"]] = A_store
    }
    
    if (has_test) {
        for (j in 1:M) {
        predictions[,, j] <- predictions[,, j] * Ysd[j] + Ymu[j]
        sigma_predictions[,, j] <- sigma_predictions[,, j] * (Ysd[j]^2)
        }
        result[["predictions"]] <- predictions
        result[["sigma_predictions"]] <- sigma_predictions
        
        # Replace your LPL aggregation with:
        LPL <- apply(LPL_draws, 2, function(v) {
        m <- max(v) # v = vector over draws for a given horizon
        m + log(mean(exp(v - m)))
        })
        names(LPL) <- paste0("t+", 1:n_ahead)
        
        # After the loop, aggregate with LSE per horizon & variable:
        LPL_univariate <- apply(PL_univariate_draws, c(2, 3), function(v) {
        m <- max(v)
        m + log(mean(exp(v - m)))
        })
        dimnames(LPL_univariate) <- list(paste0("t+", 1:n_ahead), variable_names)
        
        # numerical_normalizer <- apply(LPL_draws, 2 , max) - 700
        # LPL <- log(colMeans(exp( t(t(LPL_draws) - numerical_normalizer)))) + numerical_normalizer
        # names(LPL) <- paste0("t+", 1:n_ahead)
        #
        result$LPL <- LPL
        result$LPL_draws <- LPL_draws
        result$LPL_univariate <- LPL_univariate
    }
    
    class(result) <- "varbart"
    
    rm(forest_dataset_train_ls)
    rm(forest_model_config_mean_ls)
    rm(active_forest_mean_ls)
    rm(forest_samples_mean_ls)
    rm(global_model_config_ls)
    rm(forest_model_mean_ls)
    rm(rng)
    
    return(result)

}