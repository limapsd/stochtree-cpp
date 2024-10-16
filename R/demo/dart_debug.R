################## Simulation #################
## Main Goal: Debugg the Dirichilet Sampler ###
###############################################

library(stochtree)

# DGP:

# yt = phi * yt-1 + e_t
# yt = x*beta + e_t
# set.seed(23287)
# Tt = 201
# yt = rep(0,Tt)
# phi = 0.8
# yt[0] = rnorm(1)
# for(t in 2:Tt){
#   yt[t] = phi*yt[t-1] + rnorm(1)
# }
# yt = yt[2:Tt]


# Creating 10 noise variables

# Yt_lags = embed(yt, p+1)
# X_noise = matrix(rnorm( (Tt -2) * M), nrow = (Tt-2))
# X_test  = as.data.frame(matrix(rnorm( (Tt -2) * (M+1)), nrow = (Tt-2)))
# Yt = Yt_lags[,1]
# Xt = cbind(Yt_lags[,2],X_noise)
# Xt = as.data.frame(Xt)
#
#
#
# par(mfrow = c(1,2))
# plot(yt, main = "Yt DGP", type = "l")
# acf(yt)
M = 9
cycles = 10
probs_BART = matrix(NA, cycles, M+2)
probs_stochtree = matrix(NA, cycles, M+2)

for( i in 1:cycles){

# Tt = 200
# phi = 1.5
# xt = rnorm(Tt)
# yt = phi*xt + rnorm(Tt)
# 
# 
# 
# X_noise = matrix(rnorm( Tt * M), nrow = Tt)
# 
# Yt = yt
# 
# Xt = cbind(xt,X_noise)
# Xt = as.data.frame(Xt)

  # DGP:
# yt = phi * yt-1 + e_t
# yt = x*beta + e_t


Tt = 201
yt = rep(0,Tt)
phi = 0.8
yt[0] = rnorm(1)

for(t in 2:Tt){
  yt[t] = phi*yt[t-1] + rnorm(1)
}
yt = yt[2:Tt]
#
# Creating 10 noise variables

M = 10
p = 1

Yt_lags = embed(yt, p+1)
X_noise = matrix(rnorm( (Tt -2) * M), nrow = (Tt-2))
X_test  = as.data.frame(matrix(rnorm( (Tt -2) * (M+1)), nrow = (Tt-2)))
Yt = Yt_lags[,1]
Xt = cbind(Yt_lags[,2],X_noise)
Xt = as.data.frame(Xt)



################ Fitting the Model  ###################

num_gfr <- 0
num_burnin <- 15000
num_mcmc <- 5000
num_samples <- num_gfr + num_burnin + num_mcmc
n_trees = 200


stochtree_fit <- stochtree::bart(
  X_train = Xt, y_train = Yt,
  num_trees = n_trees, num_gfr = num_gfr, num_burnin = num_burnin,
  num_mcmc = num_mcmc, sample_sigma = T, sample_tau = F, Sparse = T, Theta_Update = T)

idx_mcmc = (num_burnin +1):num_samples
idx_dart = 1:(num_burnin/2)

variable_count_splits = stochtree_fit$model_params$variable_count_splits[idx_mcmc,]
alpha = stochtree_fit$model_params$theta_used[-idx_dart]
lpvs  = stochtree_fit$model_params$lpvs[-idx_dart,]
probs = colMeans(variable_count_splits >0)
probs_stochtree[i,] = probs

# par(mfrow = c(1,3))
# 
# plot(probs, main = "Sparse Prior", ylab = "PIP", xlab = "Predictors", ylim = c(0,1))
# plot(alpha, main = "Alpha Sampler")
# plot( alpha/(alpha + ncol(Xt) ), main = "Rho proposed", ylab = "rho")

BART_fit_dart = BART::wbart(Xt, Yt, ntree = n_trees, nskip = num_burnin, ndpost = num_mcmc, sparse = T, a=.25, rho = sqrt(ncol(Xt)))
variable_counts = BART_fit_dart$varcount
post_probs      = colMeans(variable_counts > 0)
probs_BART[i,] = post_probs


####################################################
########### SoftBart Implementation ################
####################################################

# library(SoftBart)
# SoftBART = SoftBart::softbart(Xt, Yt, X_test, opts = Opts(num_burn = num_burnin, num_save = num_mcmc, update_tau = TRUE),
#                                                         hypers = Hypers(X = Xt, Y = Yt,num_tree = n_trees))
# variable_counts = SoftBART$var_counts
# post_probs      = colMeans(variable_counts > 0)
# par(mfrow = c(1,2))
# plot(post_probs, xlab = "variables", ylab ="Posterior Probability", main = paste0("m = ",n_trees), xlim = c(1,M))
# plot(SoftBART$alpha)
# 

}

par(mfrow = c(3,2))

for(i in 1:cycles){
  plot(probs_BART[i,], main = paste0("BART- iteration  ",i))
  plot(probs_stochtree[i,], main = paste0("Stochtree- iteration  ",i))
}


