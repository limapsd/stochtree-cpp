#------------------- Aux Functions ---------------------#
upper_bound <- function(m) {
  as.integer(floor((m - 1) / 2))
}


logdmvnorm_fsv_woodbury <- function(y, mu, Lambda, h_idio, h_fac) {
  # y, mu: length M (standardized)
  # Lambda: M x Q
  # h_idio: length M, log idiosyncratic variances (HT)
  # h_fac:  length Q, log factor variances (OT)
  
  M <- length(y)
  Q <- length(h_fac)
  
  # Variances in factor model: Sigma = D + Lambda diag(f) Lambda'
  D     <- exp(h_idio)    # idiosyncratic variances
  fvar  <- exp(h_fac)     # factor variances
  
  Dinv  <- 1 / D
  Finv  <- 1 / fvar
  
  res   <- y - mu         # residual (standardized)
  L     <- Lambda
  
  ## Build K = F^{-1} + L' D^{-1} L  (Q x Q)
  # L' D^{-1} L = t(L) %*% (Dinv * L), where Dinv multiplies rows of L
  LdL   <- crossprod(L, L * Dinv)   # t(L) %*% (Dinv * L)
  K     <- diag(Finv, Q) + LdL
  
  # Cholesky of K for both logdet and solving
  Kchol <- chol(K)
  logdetK <- 2 * sum(log(diag(Kchol)))
  
  # log|Sigma| = log|D| + log|F| + log|K|
  logdetSigma <- sum(log(D)) + sum(log(fvar)) + logdetK
  
  ## Quadratic form: res' Sigma^{-1} res
  # Sigma^{-1} = D^{-1} - D^{-1} L K^{-1} L' D^{-1}
  u <- Dinv * res                     # D^{-1} res  (M x 1)
  c_vec <- crossprod(L, u)            # L' D^{-1} res  (Q x 1)
  # solve K w = c_vec
  w <- backsolve(Kchol, forwardsolve(t(Kchol), c_vec))
  quad <- sum(res * u) - sum(c_vec * w)
  
  # log density
  val <- -0.5 * (M * log(2 * pi) + logdetSigma + quad)
  as.numeric(val)
}

preprocessPredictionMatrix <- function(input_matrix, metadata) {
  # Input checks
  if (!is.matrix(input_matrix)) {
    stop("covariates provided must be a matrix")
  }
  if (!(ncol(input_matrix) == metadata$num_numeric_vars)) {
    stop(
      "Prediction set covariates have inconsistent dimension from train set covariates"
    )
  }
  
  return(input_matrix)
}
#
preprocessTrainMatrix <- function(input_matrix) {
  # Input checks
  if (!is.matrix(input_matrix)) {
    stop("covariates provided must be a matrix")
  }
  
  # Unpack metadata (assuming all variables are numeric)
  names(input_matrix) <- paste0("x", 1:ncol(input_matrix))
  df_vars <- names(input_matrix)
  num_ordered_cat_vars <- 0
  num_unordered_cat_vars <- 0
  num_numeric_vars <- ncol(input_matrix)
  numeric_vars <- names(input_matrix)
  feature_types <- rep(0, ncol(input_matrix))
  
  # Unpack data
  X <- input_matrix
  
  # Aggregate results into a list
  metadata <- list(
    feature_types = feature_types,
    num_ordered_cat_vars = num_ordered_cat_vars,
    num_unordered_cat_vars = num_unordered_cat_vars,
    num_numeric_vars = num_numeric_vars,
    numeric_vars = numeric_vars,
    original_var_indices = 1:num_numeric_vars
  )
  output <- list(
    data = X,
    metadata = metadata
  )
  
  return(output)
}

preprocessParams <- function(default_params, user_params = NULL) {
  # Override defaults from general_params
  if (!is.null(user_params)) {
    for (key in names(user_params)) {
      if (key %in% names(default_params)) {
        val <- user_params[[key]]
        if (!is.null(val)) default_params[[key]] <- val
      }
    }
  }
  
  # Return result
  return(default_params)
}

preprocessTrainData <- function(input_data) {
  # Input checks
  if ((!is.matrix(input_data)) && (!is.data.frame(input_data))) {
    stop("Covariates provided must be a dataframe or matrix")
  }
  
  # Routing the correct preprocessing function
  if (is.matrix(input_data)) {
    output <- preprocessTrainMatrix(input_data)
  } else {
    output <- preprocessTrainDataFrame(input_data)
  }
  
  return(output)
}

sample_factors <- function(resid, L, H, O) {
  q <- ncol(L)
  TT <- nrow(resid)
  F_draws <- matrix(NA_real_, TT, q)
  
  for (tt in 1:TT) {
    w <- exp(-0.5 * H[tt, ]) 
    X <- L * w 
    y <- resid[tt, ] * w 
    
    S <- crossprod(X) 
    if (q == 1L) {
      S[1, 1] <- S[1, 1] + 1.0 / exp(O[tt])
      b <- crossprod(X, y) 
    } else {
      S <- S + diag(1.0 / exp(O[tt, ]), q, q)
      b <- crossprod(X, y) 
    }
    
    R <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-6, q)))
    
    # Mean: Solve R'R mu = b
    mu <- backsolve(R, forwardsolve(t(R), b))
    
    # Sample: z ~ N(0, I), we need x ~ N(0, (R'R)^-1)
    # x = R^-1 z satisfies this.
    z <- rnorm(q)
    jitter <- backsolve(R, z) 
    
    F_draws[tt, ] <- as.numeric(mu + jitter)
  }
  F_draws
}

# resid: T x m residuals
# Ft   : T x q factors
# Ht   : T x m idiosyncratic log-vols
# V0   : either length-q prior variances for loadings, or m x q (row-specific)

sample_loadings <- function(resid, Ft, Ht, V0) {
  m <- ncol(resid)
  q <- ncol(Ft)
  Lambda_draws <- matrix(NA_real_, m, q)
  
  for (j in 1:m) {
    w <- exp(-0.5 * Ht[, j])
    X <- Ft * w 
    y <- resid[, j] * w 
    
    # V0 can be matrix (row specific) or vector
    V0_j <- if (is.matrix(V0)) V0[j, ] else V0
    
    S <- crossprod(X)
    diag(S) <- diag(S) + 1 / V0_j 
    
    R <- tryCatch(chol(S), error = function(e) chol(S + diag(1e-6, q)))
    
    b <- crossprod(X, y)
    mu <- backsolve(R, forwardsolve(t(R), b))
    
    z <- rnorm(q)
    draw <- mu + backsolve(R, z) 
    
    Lambda_draws[j, ] <- as.numeric(draw)
  }
  Lambda_draws
}


sample_hs <- function(bdraw, lambda.hs, nu.hs, tau.hs, zeta.hs) {
  k <- length(bdraw)
  if (is.na(tau.hs)) {
    tau.hs <- 1
  } else {
    tau.hs <- invgamma::rinvgamma(
      1,
      shape = (k + 1) / 2,
      rate = 1 / zeta.hs + sum(bdraw^2 / lambda.hs) / 2
    )
  }
  
  lambda.hs <- invgamma::rinvgamma(
    k,
    shape = 1.,
    rate = 1 / nu.hs + bdraw^2 / (2 * tau.hs)
  )
  
  nu.hs <- invgamma::rinvgamma(k, shape = 1., rate = 1 + 1 / lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1, shape = 1., rate = 1 + 1 / tau.hs)
  
  ret <- list(
    "psi" = (lambda.hs * tau.hs),
    "lambda" = lambda.hs,
    "tau" = tau.hs,
    "nu" = nu.hs,
    "zeta" = zeta.hs
  )
  return(ret)
}


## update one column j of Lambda
## facloads_j : length-m vector Lambda[, j]
## tau2_j     : length-m vector tau_{ij}^2 for column j
## aj, cj, dj : scalars for this column
sample_ng_col <- function(facloads_j, tau2_j, aj, cj, dj) {
  m <- length(facloads_j)

  ## 1) Sample lambda_j | tau^2
  lambda_j <- rgamma(
    n     = 1L,
    shape = cj + aj * m,
    rate  = dj + 0.5 * aj * sum(tau2_j)   # tau2_j *is* tau_{ij}^2
  )

  ## 2) Sample tau_{ij}^2 | lambda_j, Lambda_{ij}
  ## GIGrvg::rgig(n, lambda, chi, psi)
  tau2_j_new <- GIGrvg::rgig(
    n      = m,
    lambda = aj - 0.5,
    chi    = facloads_j^2,
    psi    = aj * lambda_j
  )

  list(lambda_j = lambda_j, tau2_j = tau2_j_new)
}