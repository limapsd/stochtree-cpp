
# -----------------------------------------------------------------------------------------------
# function to draw the factor loadings (basic linear regression)

get_facload <- function(yy,xx,l_sd){
  
  V_prinv <- diag(NCOL(xx))/l_sd
  V_lambda <- solve(crossprod(xx) + V_prinv)
  lambda_mean <- V_lambda %*% (crossprod(xx,yy))
  
  lambda_draw <- lambda_mean + t(chol(V_lambda)) %*% rnorm(NCOL(xx))
  return(lambda_draw)
}

# factor loadings draw
#' Title
#'
#' @param eps 
#' @param fac 
#' @param S 
#' @param pr 
#' @param m 
#' @param q 
#' @param id.fac 
#'
#' @return
#' @export
#'
#' @examples
get_Lambda <- function(eps,fac,S,pr,m,q,id.fac){
  
  L <- matrix(0,m,q)
  if(id.fac){
    for(jj in 1:m){
      if (jj<=q){
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj]-fac[,jj])/normalizer
        xx0 <- fac[,1:(jj-1),drop=FALSE]/normalizer
        if (jj>1){
          l_sd <- pr[jj,1:(jj-1)]
          lambda0 <- get_facload(yy0,xx0,l_sd=l_sd)
        }else{
          lambda0 <- 1
        }
        
        if (jj>1){
          L[jj,1:(jj-1)] <- lambda0
          L[jj,jj] <- 1
        }else if (jj==1){
          L[jj,jj] <- 1
        }
      }else{
        normalizer <- exp(0.5*S[,jj])
        yy0 <- (eps[,jj])/normalizer
        xx0 <- fac[,,drop=FALSE]/normalizer
        l_sd <- pr[jj,]
        lambda0 <- get_facload(yy0,xx0,l_sd=l_sd)
        L[jj,] <- lambda0
      }
    }
  }else{
    for(jj in 1:m){
      normalizer <- exp(0.5*S[,jj])
      yy0 <- (eps[,jj])/normalizer
      xx0 <- fac[,,drop=FALSE]/normalizer
      l_sd <- pr[jj,]
      lambda0 <- get_facload(yy0,xx0,l_sd=l_sd)
      L[jj,] <- lambda0
    }
  }
  return(L)
}

# sample the latent factors
#' Title
#'
#' @param e 
#' @param S 
#' @param H 
#' @param L 
#' @param q 
#' @param t 
#'
#' @return
#' @export
#'
#' @examples
get_factors <- function(e,S,H,L,q,t){
  
  F_raw <- matrix(0,t,q)
  for (tt in 1:t){
    normalizer <- exp(-S[tt,]/2)
    Lt <- L*normalizer
    yt <- e[tt,]*normalizer
    
    if (q==1) fac.varinv <-  1/H[tt] else fac.varinv <- diag(q)/H[tt]
    fac.Sigma <-  solve(crossprod(Lt)+fac.varinv)
    fac.mean <- fac.Sigma%*%crossprod(Lt,yt)
    
    F_temp <- try(fac.mean + t(chol(fac.Sigma)) %*% rnorm(q),silent=TRUE)
    if (is(F_temp,"try-error")) F_temp <- fac.mean + t(chol(fac.Sigma+diag(q)*1e-6)) %*% rnorm(q)
    F_raw[tt,] <- F_temp
  }
  return(F_raw)
}


# -----------------------------------------------------------------------------------------------
# function to draw the offdiagnoal horseshoe prior 
## Sample HS for a regression.

#' Title
#'
#' @param bdraw 
#' @param lambda.hs 
#' @param nu.hs 
#' @param tau.hs 
#' @param zeta.hs 
#'
#' @return
#' @export
#'
#' @examples
get_hs <- function(bdraw,lambda.hs,nu.hs,tau.hs,zeta.hs){
  k <- length(bdraw)
  if (is.na(tau.hs)){
    tau.hs <- 1   
  }else{
    tau.hs <- invgamma::rinvgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lambda.hs)/2) 
  }
  
  lambda.hs <- invgamma::rinvgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
  
  nu.hs <- invgamma::rinvgamma(k,shape=1,rate=1+1/lambda.hs)
  zeta.hs <- invgamma::rinvgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lambda.hs*tau.hs),"lambda"=lambda.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}
