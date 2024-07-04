library("bayesianVARs")
library("stochtree")
# source("flexBART2.R")
# source("varbartfunc.R")
# source("varbart.R")

###############################################
###############################################
############### Lag Function ##################

# lag variables
mlag <- function(X,lag){
  p <- lag
  X <- as.matrix(X)
  Traw <- nrow(X)
  N <- ncol(X)
  Xlag <- matrix(NA,Traw,p*N)
  for (ii in 1:p){
    Xlag[(p+1):Traw,(N*(ii-1)+1):(N*ii)]=X[(p+1-ii):(Traw-ii),(1:N)]
  }
  return(Xlag)
}

#############################################
#############################################

p = 10
# variables = c( "GDPC1","CPIAUCSL")

variables = c("GDPC1","INDPRO", "CPIAUCSL","FEDFUNDS","GS10")
train_data = 100 * bayesianVARs::usmacro_growth[1:231, variables]
Yraw = train_data 

Ymu  = apply(Yraw, 2, mean,na.rm=T)
Ysd  = apply(Yraw, 2, sd,na.rm=T)
Yraw_stdize = apply(Yraw, 2, function(x){(x-mean(x,na.rm=T))/sd(x,na.rm=T)})

# create design matrices Y/X

X_train = cbind(mlag(Yraw_stdize,p))[(p+1):nrow(Yraw_stdize),]
Y_train = Yraw_stdize[(p+1):nrow(Yraw_stdize), ] # leave original input matrix unchanged

varbart_model = stochtree::varbart(Y_train, X_train, Yraw, p, 5, Heteroskedacity_Model = "CSV", num_burnin = 1000, num_mcmc = 1000)
fitted_model = varbart_model$y_hat_insample
vol_model = varbart_model$H_matrix
splits = varbart_model$variable_count_splits

post_probs <- colMeans(splits > 0)
plot(post_probs[,1], ylim = c(0,1), xlab = "variables", ylab ="Posterior Probability", main = paste0("m = ",250))


par(mfrow = c(2,2))
for( i in 1:length(variables)){
  var_names = variables[i]
  bart.h.fit = fitted_model[,,i]
  plot(apply(sqrt(vol_model[,,i]),2, median), ylab = expression(sigma), main = paste0(var_names, "- Stochtree Implementation"))
  quantiles.bart.h = apply(bart.h.fit, 2 , quantile, probs = c(0.05,0.5,0.95))
  print(mean((quantiles.bart.h["50%",] - train_data[-1:-p,i])^2))
  plot(train_data[-1:-p,i],type = "l", pch=16, cex=0.75, xlab = "pred", ylab = "actual", main =paste0(var_names," - Stochtree Implementation"))

  lines(quantiles.bart.h["50%",], col = "red")
  # lines(quantiles.bart.h["5%",], lty = "dashed", col =2)
  # lines(quantiles.bart.h["95%",], lty = "dashed", col =2)
}


fsv_varbart_model = stochtree::fsv_varbart(Y_train, X_train, Yraw, p, 5, Heteroskedacity_Model = "FSV", num_burnin = 1000, num_mcmc = 1000)
fitted_model = fsv_varbart_model$y_hat_insample
vol_model = fsv_varbart_model$H_matrix

par(mfrow = c(2,2))
for( i in 1:length(variables)){
  var_names = variables[i]
  bart.h.fit = fitted_model[,,i]
  plot(apply(sqrt(vol_model[,,i]),2, median), ylab = expression(sigma), main = paste0(var_names, "- Stochtree Implementation"))
  quantiles.bart.h = apply(bart.h.fit, 2 , quantile, probs = c(0.05,0.5,0.95))
  print(mean((quantiles.bart.h["50%",] - train_data[-1:-p,i])^2))
  plot(train_data[-1:-p,i],type = "l", pch=16, cex=0.75, xlab = "pred", ylab = "actual", main =paste0(var_names," - Stochtree Implementation"))

  lines(quantiles.bart.h["50%",], col = "red")
  lines(quantiles.bart.h["5%",], lty = "dashed", col =2)
  lines(quantiles.bart.h["95%",], lty = "dashed", col =2)
  }



varbart_model = stochtree::varbart(Y_train, X_train, Yraw, p, 5, Heteroskedacity_Model = "NONE", num_burnin = 1000, num_mcmc = 1000)
fitted_model = varbart_model$y_hat_insample
vol_model = varbart_model$H_matrix

par(mfrow = c(2,2))
for( i in 1:length(variables)){
  var_names = variables[i]
  bart.h.fit = fitted_model[,,i]
  plot(apply(vol_model[,,i],1, median), ylab = expression(sigma), main = paste0(var_names, "- Stochtree Implementation"))
  quantiles.bart.h = apply(bart.h.fit, 2 , quantile, probs = c(0.05,0.5,0.95))
  print(mean((quantiles.bart.h["50%",] - train_data[-1:-p,i])^2))
  plot(train_data[-1:-p,i],type = "l", pch=16, cex=0.75, xlab = "pred", ylab = "actual", main =paste0(var_names," - Stochtree Implementation"))

  lines(quantiles.bart.h["50%",], col = "red")
  lines(quantiles.bart.h["5%",], lty = "dashed", col ="green")
  lines(quantiles.bart.h["95%",], lty = "dashed", col ="green")
}

