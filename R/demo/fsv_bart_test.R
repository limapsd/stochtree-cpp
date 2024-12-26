library("bayesianVARs")
library("stochtree")
# source("flexBART2.R")
# source("varbartfunc.R")
# source("varbart.R")


#############################################
#############################################

p = 5
variables = c( "GDPC1","CPIAUCSL")

# variables = c("GDPC1","CPIAUCSL","FEDFUNDS","GS10")

train_data = 100 * bayesianVARs::usmacro_growth[1:231, variables]
test_data  = 100 * bayesianVARs::usmacro_growth[232:247, variables]

Yraw = train_data 


fsv_dart_model = stochtree::fsv_varbart(Yraw, p, 2, bart_prior = "dart", sv = "no", num_burnin = 20000, num_mcmc = 5000)








# par(mfrow = c(1,3))
# for( i in 1:length(variables)){
#   var_names = variables[i]
#   bart.h.fit = fitted_model[,,i]
#   plot(apply(sqrt(vol_model[,,i]),2, median), ylab = expression(sigma), main = paste0(var_names, "- Stochtree Implementation"))
#   quantiles.bart.h = apply(bart.h.fit, 2 , quantile, probs = c(0.05,0.5,0.95))
#   print( sqrt(mean((quantiles.bart.h["50%",] - train_data[-1:-p,i])^2)))
#   plot(train_data[-1:-p,i],type = "l", pch=16, cex=0.75, xlab = "t", ylab = "yt", main =paste0(var_names," - Stochtree Implementation"))
#   
#   lines(quantiles.bart.h["50%",], col = "red")
#   
#   x_vals <- seq_along(quantiles.bart.h["50%", ])  # Sequence for x-axis
#   
#   polygon(
#     c(x_vals, rev(x_vals)),  # X-coordinates for the polygon
#     c(quantiles.bart.h["95%", ], rev(quantiles.bart.h["5%", ])),  # Y-coordinates
#     col = rgb(1, 0, 0, 0.5),  # Transparent red
#     border = NA
#   )
#   
#   vc_splits = variable_count_splits_matrix[,,i]
#   plot(colMeans(vc_splits)>0, ylab = "PIP", main = paste0(var_names, " Split Probs"), ylim = c(0,1))
# }
