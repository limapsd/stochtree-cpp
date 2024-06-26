% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bart.R
\name{bart}
\alias{bart}
\title{Run the BART algorithm for supervised learning.}
\usage{
bart(
  X_train,
  y_train,
  W_train = NULL,
  group_ids_train = NULL,
  rfx_basis_train = NULL,
  X_test = NULL,
  W_test = NULL,
  group_ids_test = NULL,
  rfx_basis_test = NULL,
  cutpoint_grid_size = 100,
  tau_init = NULL,
  alpha = 0.95,
  beta = 2,
  min_samples_leaf = 5,
  leaf_model = 0,
  nu = 3,
  lambda = NULL,
  a_leaf = 3,
  b_leaf = NULL,
  q = 0.9,
  sigma2_init = NULL,
  num_trees = 200,
  num_gfr = 5,
  num_burnin = 0,
  num_mcmc = 100,
  sample_sigma = T,
  sample_tau = T,
  random_seed = -1,
  keep_burnin = F,
  keep_gfr = F,
  Sparse = F,
  Theta_Update = F
)
}
\arguments{
\item{X_train}{Covariates used to split trees in the ensemble.}

\item{y_train}{Outcome to be modeled by the ensemble.}

\item{W_train}{(Optional) Bases used to define a regression model \code{y ~ W} in
each leaf of each regression tree. By default, BART assumes constant leaf node
parameters, implicitly regressing on a constant basis of ones (i.e. \code{y ~ 1}).}

\item{group_ids_train}{(Optional) Group labels used for an additive random effects model.}

\item{rfx_basis_train}{(Optional) Basis for "random-slope" regression in an additive random effects model.
If \code{group_ids_train} is provided with a regression basis, an intercept-only random effects model
will be estimated.}

\item{X_test}{(Optional) Test set of covariates used to define "out of sample" evaluation data.}

\item{W_test}{(Optional) Test set of bases used to define "out of sample" evaluation data.
While a test set is optional, the structure of any provided test set must match that
of the training set (i.e. if both X_train and W_train are provided, then a test set must
consist of X_test and W_test with the same number of columns).}

\item{group_ids_test}{(Optional) Test set group labels used for an additive random effects model.
We do not currently support (but plan to in the near future), test set evaluation for group labels
that were not in the training set.}

\item{rfx_basis_test}{(Optional) Test set basis for "random-slope" regression in additive random effects model.}

\item{cutpoint_grid_size}{Maximum size of the "grid" of potential cutpoints to consider. Default: 100.}

\item{tau_init}{Starting value of leaf node scale parameter. Calibrated internally as \code{1/num_trees} if not set here.}

\item{alpha}{Prior probability of splitting for a tree of depth 0. Tree split prior combines \code{alpha} and \code{beta} via \code{alpha*(1+node_depth)^-beta}.}

\item{beta}{Exponent that decreases split probabilities for nodes of depth > 0. Tree split prior combines \code{alpha} and \code{beta} via \code{alpha*(1+node_depth)^-beta}.}

\item{min_samples_leaf}{Minimum allowable size of a leaf, in terms of training samples. Default: 5.}

\item{nu}{Shape parameter in the \code{IG(nu, nu*lambda)} global error variance model. Default: 3.}

\item{lambda}{Component of the scale parameter in the \code{IG(nu, nu*lambda)} global error variance prior. If not specified, this is calibrated as in Sparapani et al (2021).}

\item{a_leaf}{Shape parameter in the \code{IG(a_leaf, b_leaf)} leaf node parameter variance model. Default: 3.}

\item{b_leaf}{Scale parameter in the \code{IG(a_leaf, b_leaf)} leaf node parameter variance model. Calibrated internally as 0.5/num_trees if not set here.}

\item{q}{Quantile used to calibrated \code{lambda} as in Sparapani et al (2021). Default: 0.9.}

\item{sigma2_init}{Starting value of global variance parameter. Calibrated internally as in Sparapani et al (2021) if not set here.}

\item{num_trees}{Number of trees in the ensemble. Default: 200.}

\item{num_gfr}{Number of "warm-start" iterations run using the grow-from-root algorithm (He and Hahn, 2021). Default: 5.}

\item{num_burnin}{Number of "burn-in" iterations of the MCMC sampler. Default: 0.}

\item{num_mcmc}{Number of "retained" iterations of the MCMC sampler. Default: 100.}

\item{sample_sigma}{Whether or not to update the \code{sigma^2} global error variance parameter based on \code{IG(nu, nu*lambda)}. Default: T.}

\item{sample_tau}{Whether or not to update the \code{tau} leaf scale variance parameter based on \code{IG(a_leaf, b_leaf)}. Cannot (currently) be set to true if \code{ncol(W_train)>1}. Default: T.}

\item{random_seed}{Integer parameterizing the C++ random number generator. If not specified, the C++ random number generator is seeded according to \code{std::random_device}.}

\item{keep_burnin}{Whether or not "burnin" samples should be included in cached predictions. Default FALSE. Ignored if num_mcmc = 0.}

\item{keep_gfr}{Whether or not "grow-from-root" samples should be included in cached predictions. Default TRUE. Ignored if num_mcmc = 0.}

\item{Sparse}{Whether you want to turn on the dirichilet prior.}

\item{Theta_Update}{Whether or not update the theta of the dirichilet prior.}

\item{variable_weights}{Vector of length \code{ncol(X_train)} indicating a "weight" placed on each
variable for sampling purposes. Default: \code{rep(1/ncol(X_train),ncol(X_train))}.}
}
\value{
List of sampling outputs and a wrapper around the sampled forests (which can be used for in-memory prediction on new data, or serialized to JSON on disk).
}
\description{
Run the BART algorithm for supervised learning.
}
\examples{
n <- 100
p <- 5
X <- matrix(runif(n*p), ncol = p)
f_XW <- (
    ((0 <= X[,1]) & (0.25 > X[,1])) * (-7.5) + 
    ((0.25 <= X[,1]) & (0.5 > X[,1])) * (-2.5) + 
    ((0.5 <= X[,1]) & (0.75 > X[,1])) * (2.5) + 
    ((0.75 <= X[,1]) & (1 > X[,1])) * (7.5)
)
noise_sd <- 1
y <- f_XW + rnorm(n, 0, noise_sd)
test_set_pct <- 0.2
n_test <- round(test_set_pct*n)
n_train <- n - n_test
test_inds <- sort(sample(1:n, n_test, replace = F))
train_inds <- (1:n)[!((1:n) \%in\% test_inds)]
X_test <- X[test_inds,]
X_train <- X[train_inds,]
y_test <- y[test_inds]
y_train <- y[train_inds]
bart_model <- bart(X_train = X_train, y_train = y_train, X_test = X_test)
# plot(rowMeans(bart_model$y_hat_test), y_test, xlab = "predicted", ylab = "actual")
# abline(0,1,col="red",lty=3,lwd=3)
}
