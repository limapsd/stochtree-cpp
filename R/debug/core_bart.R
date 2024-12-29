######### This is the core implementation of a BART algo ###############

# Generate the data
n <- 500
p_X <- 10
p_W <- 1
X <- matrix(runif(n*p_X), ncol = p_X)
f_XW <- (
  ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) +
    ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) +
    ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) +
    ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
)
y <- f_XW + rnorm(n, 0, 1)

# Standardize outcome
y_bar <- mean(y)
y_std <- sd(y)
resid <- (y-y_bar)/y_std
####################################################################

alpha <- 0.95 # standard
beta <- 2.0  # standard
min_samples_leaf <- 5
max_depth <- 10
n_trees <- 200
cutpoint_grid_size = 100
global_variance_init = 1.
tau_init = 0.5
leaf_prior_scale = matrix(c(tau_init), ncol = 1)
nu <- 4
lambda <- 0.5
a_leaf <- 2.
b_leaf <- 0.5

feature_types <- as.integer(rep(0, p_X)) # 0 = numeric
var_weights <- rep(1/p_X, p_X)
outcome_model_type <- 0


# Creating External Pointers:

outcome <- createOutcome(resid)
rng <- createRNG()
forest_dataset <- createForestDataset(X)

forest_model <- createForestModel(forest_dataset, feature_types,
                                  n_trees, n, alpha, beta,
                                  min_samples_leaf, max_depth)
forest_samples <-createForestContainer(num_trees = n_trees, is_leaf_constant = T)
active_forest  <-createForest(num_trees = n_trees, is_leaf_constant = T)

num_gfr <- 0
num_mcmc   <- 10000
num_burnin <- 10000
num_samples <- num_burnin + num_mcmc
keep_burnin <- FALSE

global_var_samples <- c(global_variance_init, rep(0, num_samples))
leaf_scale_samples <- c(tau_init, rep(0, num_samples))

for (i in 1:num_samples) {
  is_mcmc <- i > num_burnin
  if (is_mcmc) {
    mcmc_counter <- i - (num_gfr + num_burnin)
    if (mcmc_counter %% 1 == 0) keep_sample <- T
    else keep_sample <- F
  } else {
    if (keep_burnin) keep_sample <- T
    else keep_sample <- F
  }
  
  # Sample forest
  forest_model$sample_one_iteration(
    forest_dataset, outcome, forest_samples, active_forest, rng, feature_types,
    outcome_model_type, leaf_prior_scale, var_weights,
    1, 1, global_var_samples[i], cutpoint_grid_size, keep_forest = keep_sample, gfr = F
  )
  
  # Sample global variance parameter
  global_var_samples[i+1] <- sample_sigma2_one_iteration(
    outcome, forest_dataset, rng, nu, lambda
  )
  
  # Sample leaf node variance parameter and update `leaf_prior_scale`
  leaf_scale_samples[i+1] <- sample_tau_one_iteration(
    active_forest, rng, a_leaf, b_leaf
  )
  leaf_prior_scale[1,1] <- leaf_scale_samples[i+1]
}
active_forest$get_forest_split_counts(ncol(X))
# Forest predictions

preds <- forest_samples$predict(forest_dataset)*y_std + y_bar

# Global error variance
sigma_samples <- sqrt(global_var_samples)*y_std

plot(sigma_samples[(num_burnin+1):num_samples], ylab="sigma")
abline(h = 1, lty = 2, col = 2)

plot(rowMeans(preds[,1:num_mcmc]), y, pch=16,
     cex=0.75, xlab = "pred", ylab = "actual")
abline(0,1,col="red",lty=2,lwd=2.5)

rm(forest_model)
rm(forest_dataset)
rm(outcome)
rm(rng)




