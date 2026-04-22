#' Generalized Impulse Response Function for VARBART
#'
#' Computes Koop-Pesaran-Potter (1996) generalized impulse responses for a
#' fitted VARBART model. Supports the homoskedastic case
#' (variance_prior = "const") only. The model must have been fit with
#' return_forests = TRUE.
#'
#' Convention: the shock occurs at the impact period, h = 1 in the returned
#' arrays (i.e. the first horizon is the impact response). The shock REPLACES
#' the residual draw of variable shock_var at h = 1 (KPP-style), keeping all
#' other contemporaneous and future shocks identical between baseline and
#' shocked paths (common-shock variance reduction).
#'
#' @param model      A varbart object fit with variance_prior = "const",
#'                   mean_prior = "cgm", and return_forests = TRUE.
#' @param history    A p x M matrix in ORIGINAL units, with the most-recent
#'                   lag (y_{t-1}) in row 1, y_{t-2} in row 2, etc.
#' @param shock_var  Index in 1..M of the variable to shock.
#' @param shock_size Size of the shock in ORIGINAL units (added to the
#'                   reduced-form residual of variable shock_var).
#' @param horizon    Number of horizons to compute.
#' @param R          Monte Carlo paths per posterior draw (variance reduction).
#' @param mc_draws   Optional vector of MCMC draw indices to use; defaults to all.
#' @param verbose    Show progress bar.
#' @returns A list with:
#'   - girf_draws: array (num_draws x horizon x M) of GIRFs in original scale
#'   - median:     horizon x M pointwise median
#'   - quantiles:  4 x horizon x M pointwise quantiles (5%, 16%, 84%, 95%)
#' @export
#'
#' @examples
girf <- function(model,
                 history,
                 shock_var,
                 shock_size,
                 horizon = 12,
                 R = 200,
                 mc_draws = NULL,
                 verbose = TRUE) {

  ## ---- Validation ---------------------------------------------------------
  if (is.null(model$forest_samples_mean_ls)) {
    stop("Model was not fit with return_forests = TRUE.")
  }
  vp <- model$model_params$variance_prior
  if (!(vp %in% c("const", "csv"))) {
    stop("This GIRF implementation supports variance_prior in {'const','csv'}.")
  }

  M   <- model$M
  p   <- model$lags
  K   <- model$K
  Ymu <- model$model_params$outcome_mean
  Ysd <- model$model_params$outcome_scale

  if (!is.matrix(history) || nrow(history) != p || ncol(history) != M) {
    stop(sprintf("history must be a %d x %d matrix (most-recent lag row 1).", p, M))
  }
  if (!(shock_var %in% seq_len(M))) {
    stop(sprintf("shock_var must be in 1..%d", M))
  }

  if (vp == "const") {
    if (is.null(model$sigma_store)) {
      stop("model$sigma_store missing. Use the patched varbart() that stores it.")
    }
    sigma2_orig <- model$sigma_store[, 1, , drop = TRUE]
    if (is.null(dim(sigma2_orig))) sigma2_orig <- matrix(sigma2_orig, nrow = 1)
    sigma2_std  <- sweep(sigma2_orig, 2, Ysd^2, "/")
  } else if (vp == "csv") {
    if (is.null(model$A_store) || is.null(model$Ht_store) ||
        is.null(model$sv_params_mcmc)) {
      stop("model$A_store / Ht_store / sv_params_mcmc missing for CSV.")
    }
  }

  num_draws_total <- dim(model$y_hat_train)[1]
  if (is.null(mc_draws)) mc_draws <- seq_len(num_draws_total)
  num_draws <- length(mc_draws)

  history_std <- sweep(sweep(history, 2, Ymu, "-"), 2, Ysd, "/")
  x0_std      <- as.numeric(t(history_std))
  shock_std   <- shock_size / Ysd[shock_var]

  girf_draws <- array(0, dim = c(num_draws, horizon, M))

  if (verbose) {
    cat(sprintf(
      "GIRF [%s]: %d draws x %d paths x %d horizons x %d eqns\n",
      vp, num_draws, R, horizon, M
    ))
    pb <- txtProgressBar(min = 0, max = num_draws, style = 3)
  }

  for (idx_m in seq_len(num_draws)) {
    m <- mc_draws[idx_m]

    x_state_base    <- matrix(rep(x0_std, each = R), nrow = R)
    x_state_shocked <- x_state_base

    if (vp == "const") {
      sd_m <- sqrt(as.numeric(sigma2_std[m, ]))
    } else {
      A_m    <- model$A_store[m, , ]
      Ainv_m <- forwardsolve(A_m, diag(M))

      TT_store <- dim(model$Ht_store)[2]
      h_init   <- log(as.numeric(model$Ht_store[m, TT_store, ]))

      mu_m    <- as.numeric(model$sv_params_mcmc[m, 1, ])
      phi_m   <- as.numeric(model$sv_params_mcmc[m, 2, ])
      sigma_m <- as.numeric(model$sv_params_mcmc[m, 3, ])

      h_state <- matrix(h_init, R, M, byrow = TRUE)
    }

    diff_path <- array(0, dim = c(R, horizon, M))

    for (h in seq_len(horizon)) {

      if (vp == "const") {
        eps_base    <- matrix(rnorm(R * M), R, M) *
                       matrix(sd_m, R, M, byrow = TRUE)
        eps_shocked <- eps_base
        if (h == 1L) eps_shocked[, shock_var] <- shock_std
        eta_base    <- eps_base
        eta_shocked <- eps_shocked
      } else {
        xi <- matrix(rnorm(R * M), R, M)
        h_new <- matrix(mu_m,    R, M, byrow = TRUE) +
                 matrix(phi_m,   R, M, byrow = TRUE) *
                   (h_state - matrix(mu_m, R, M, byrow = TRUE)) +
                 matrix(sigma_m, R, M, byrow = TRUE) * xi

        sd_paths    <- sqrt(exp(h_new))
        eps_base    <- matrix(rnorm(R * M), R, M) * sd_paths
        eps_shocked <- eps_base
        if (h == 1L) eps_shocked[, shock_var] <- shock_std

        eta_base    <- eps_base    %*% t(Ainv_m)
        eta_shocked <- eps_shocked %*% t(Ainv_m)

        h_state <- h_new
      }

      x_stack  <- rbind(x_state_base, x_state_shocked)
      ds_stack <- stochtree::createForestDataset(x_stack)

      g_base    <- matrix(0, R, M)
      g_shocked <- matrix(0, R, M)
      for (mm in seq_len(M)) {
        pred <- .predict_from_draw(
          model$forest_samples_mean_ls[[mm]], ds_stack, m
        )
        g_base[, mm]    <- pred[1:R]
        g_shocked[, mm] <- pred[(R + 1):(2 * R)]
      }

      y_h_base    <- g_base    + eta_base
      y_h_shocked <- g_shocked + eta_shocked

      diff_path[, h, ] <- y_h_shocked - y_h_base

      if (p > 1L) {
        x_state_base <- cbind(
          y_h_base,
          x_state_base[, 1:((p - 1) * M), drop = FALSE]
        )
        x_state_shocked <- cbind(
          y_h_shocked,
          x_state_shocked[, 1:((p - 1) * M), drop = FALSE]
        )
      } else {
        x_state_base    <- y_h_base
        x_state_shocked <- y_h_shocked
      }
    }

    girf_draws[idx_m, , ] <- apply(diff_path, c(2, 3), mean)

    if (verbose) setTxtProgressBar(pb, idx_m)
  }
  if (verbose) close(pb)

  girf_orig <- sweep(girf_draws, 3, Ysd, "*")

  list(
    girf_draws = girf_orig,
    median     = apply(girf_orig, c(2, 3), median),
    quantiles  = apply(girf_orig, c(2, 3), quantile,
                       probs = c(0.05, 0.16, 0.84, 0.95))
  )
}

## --- internal: predict from a single saved forest ---------------------------
.girf_env <- new.env(parent = emptyenv())
.girf_env$slow_warned <- FALSE

.predict_from_draw <- function(forest_samples, dataset, m) {
  fnum <- as.integer(m - 1L)
  res <- tryCatch(
    forest_samples$predict_raw_single_forest(dataset, fnum),
    error = function(e) NULL
  )
  if (!is.null(res)) return(as.numeric(res))

  if (!.girf_env$slow_warned) {
    warning(
      "predict_raw_single_forest not available; falling back to bulk predict ",
      "+ column extraction. GIRF will be ~num_saved times slower than needed.",
      call. = FALSE, immediate. = TRUE
    )
    .girf_env$slow_warned <- TRUE
  }
  full <- forest_samples$predict(dataset)
  as.numeric(full[, m])
}
