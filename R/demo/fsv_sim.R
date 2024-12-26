##################################
#---Simulating an VAR(2)-FSV-----#
##################################
library(MASS)
library(stochvol)
library(factorstochvol)

Tt = 400
n = 10
k = 2
p = 2

Phi_1 = matrix(0,nrow = n, ncol = n)
Phi_2 = matrix(0, n, n)
mu = matrix(runif(n, -10,10), ncol= 1)


Yt = matrix(0,nrow = Tt, ncol = n)
Yt[1:p,] =rnorm(p)


diag(Phi_1) = runif(n,-0.2,0.4)
# diag(Phi_1) = 2
diag(Phi_2) = rnorm(n, mean = 0, sd = 0.05)

for(i in 1:n){
  for(j in 1:n){
    if(j !=i){
      Phi_1[i,j] = runif(1, -0.2,0.2)
      Phi_2[i,j] = rnorm(1, 0, 0.05)
    }
  }
}

idipara = matrix(0,nrow = n , ncol = 3) #mu, phi, sigma
idipara[,1] = rep(-1,n)
idipara[,2] = rep(0.98,n)
idipara[,3] = rep(0.317,n)

facpara = matrix(0, nrow =k, ncol =2) #phi, sigma
facpara[,1] = rep(0.98,k)
facpara[,2] = rep(0.317,k)

et = factorstochvol::fsvsim(n = Tt, series = n, factors = k, idipara = idipara, facpara = facpara)
e_t = et$y

for(t in (p+1):Tt){
  Yt[t,] = mu + Phi_1%*%Yt[t-1,] + Phi_2%*%Yt[t-2,] + e_t[t,]
}

ts.plot(Yt, col = 1:n)


vbart_fsv = stochtree::fsv_varbart(Yt, 2,2,sv="SV", num_burnin = 20000, num_mcmc = 10000)
yt_2 = colMeans(vbart_fsv$y_hat_insample[,,2])
plot(Yt[-1:-p,2],type="l")
lines(yt_2, col =2)


prior_phi <- specify_prior_phi(data = Yt,
                               lags = p,
                               prior = "normal")

prior_sigma <- specify_prior_sigma(data = Yt,
                                   type = "cholesky",
                                   cholesky_heteroscedastic = T,
                                   cholesky_U_prior = "HS")

mod <- bvar(Yt, lags = p, draws = 1000, burnin = 1000,
            prior_phi = prior_phi, prior_sigma = prior_sigma,
            sv_keep = "all")  

fitted_mod = fitted(mod)
fitted_model = apply(fitted_mod$fitted, 1:2, mean)
vol_model = exp(apply(mod$logvar,1:2,mean))

plot(Yt[-1:-p,2],type="l")
lines(fitted_model[,2], col =2)
lines(yt_2, col = "green")

par(mfrow = c(1,2))
plot(Yt[-1:-p,2] ~ yt_2, ylab = "Yt", xlab = "Yt_hat", xlim = c(min(fitted_model[,2]), max(fitted_model[,2])))    
abline(0,1,col="red",lty=2,lwd=2.5)

plot(Yt[-1:-p,2] ~ fitted_model[,2], ylab = "Yt", xlab = "Yt_hat")    
abline(0,1,col="red",lty=2,lwd=2.5)

cat("VAR:",sqrt(mean((Yt[-1:-p,2] - fitted_model[,2])^2)))
cat("BART", sqrt(mean((Yt[-1:-p,2] - yt_2)^2) ) )
