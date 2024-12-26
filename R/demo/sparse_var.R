M = 10
Tt = 201
cycles = 15
probs_stochtree = matrix(NA, cycles, M)

for( i in 1:cycles){
Yt = matrix(0,nrow = Tt,ncol =M)
Phi = matrix(0,nrow = M, ncol= M)
diag(Phi) = 0.8
 
Yt[1,] = rnorm(M)

for(t in 2:Tt){
Yt[t, ] = Yt[t-1, ]%*%t(Phi) + rnorm(M)
}


Yt_lags = embed(Yt,2)

Yt = Yt_lags[,1,drop = F]
Xt = Yt_lags[,(M+1):ncol(Yt_lags), drop=F]
Xt = as.data.frame(Xt)

num_gfr <- 0
num_burnin <- 20000
num_mcmc <- 5000
num_samples <- num_gfr + num_burnin + num_mcmc
n_trees = 200


stochtree_fit <- stochtree::bart(
  X_train = Xt, y_train = Yt,
  num_trees = n_trees, num_gfr = num_gfr, num_burnin = num_burnin,
  num_mcmc = num_mcmc, sample_sigma = T, sample_tau = T, Sparse = T, Theta_Update = T)

idx_mcmc = (num_burnin +1):num_samples
idx_dart = 1:(num_burnin/2)

variable_count_splits = stochtree_fit$model_params$variable_count_splits[idx_mcmc,]
alpha = stochtree_fit$model_params$theta_used[-idx_dart]
lpvs  = stochtree_fit$model_params$lpvs[-idx_dart,]
probs = colMeans(variable_count_splits >0)
probs_stochtree[i,] = probs

par(mfrow = c(1,3))

plot(probs, main = "Sparse Prior", ylab = "PIP", xlab = "Predictors", ylim = c(0,1))
plot(alpha, main = "Alpha Sampler")
plot( alpha/(alpha + ncol(Xt) ), main = "Rho proposed", ylab = "rho")
}

par(mfrow = c(3,2))

for(i in 1:cycles){
  # plot(probs_BART[i,], main = paste0("BART- iteration  ",i))
  plot(probs_stochtree[i,], main = paste0("Stochtree- iteration  ",i), ylim = c(0,1))
}
