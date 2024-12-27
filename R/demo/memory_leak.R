library(stochtree)

# Generate the data
n <- 500
p_x <- 10
snr <- 3
X <- matrix(runif(n*p_x), ncol = p_x)
f_XW <- (
  ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) + 
    ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) + 
    ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) + 
    ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
)
noise_sd <- sd(f_XW) / snr
y <- f_XW + rnorm(n, 0, 1)*noise_sd

# Split data into test and train sets
test_set_pct <- 0.2
n_test <- round(test_set_pct*n)
n_train <- n - n_test
test_inds <- sort(sample(1:n, n_test, replace = FALSE))
train_inds <- (1:n)[!((1:n) %in% test_inds)]
X_test <- as.data.frame(X[test_inds,])
X_train <- as.data.frame(X[train_inds,])
W_test <- NULL
W_train <- NULL
y_test <- y[test_inds]
y_train <- y[train_inds]


num_gfr <- 0
num_burnin <- 30000
num_mcmc <- 10000
num_samples <- num_gfr + num_burnin + num_mcmc
n_trees <- 200

stochtree_fit <- stochtree::bart(
  X_train = X_train, y_train = y_train,
  num_trees = n_trees, num_gfr = num_gfr, num_burnin = num_burnin,
  num_mcmc = num_mcmc, sample_sigma = T, sample_tau = F, Sparse = F, Theta_Update = F)
