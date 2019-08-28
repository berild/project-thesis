library(INLA)
library(tidyverse)

draw.prop.beta <- function(b,sigma = 0.75){
  rnorm(ncol(b),mean = as.numeric(b), sd = sigma)
}

prob.prop.beta <- function(b1,b2, sigma = 0.75, log = TRUE){
  sum(dnorm(as.numeric(b1),mean = as.numeric(b2), sd = sigma, log = log))
}

prior.beta <- function(b, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(as.numeric(b), mean = 0, sd= sigma, log = log))
}


fit.inla <- function(data,b){
  data$oset = b$beta_1 * data$x1 +b$beta_2*data$x2
  res = inla(y ~1+offset(oset), data = data)
  return(list(mlik = res$mlik[[1]], alfa = res$marginals.fixed[[1]], tau = res$marginals.hyperpar[[1]]))
}


sample.data.bi.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2
  df = data.frame(y = y, x1 = x1, x2 = x2)
}


bi.linreg.mcmc.w.inla <- function(df){
  N = 10000
  burnin = 500
  beta = data.frame(beta_1 = rep(NA,N),beta_2 = rep(NA,N),is_burnin = c(rep(T,burnin),rep(F,N-burnin)))
  beta[1,-3] = c(0,0)
  mlik.y1 = fit.inla(df,beta[1,])
  first = T
  browser()
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2, N)){
    setTxtProgressBar(pb, i)
    beta[i,-3] = draw.prop.beta(beta[i-1,])
    mlik.y2 = fit.inla(df,beta[i,])
    acc.prob = c(acc.prob,mlik.y2$mlik + prior.beta(beta[i,]) +
      prob.prop.beta(beta[i-1,],beta[i,]) - mlik.y1$mlik -
      prior.beta(beta[i-1,]) - prob.prop.beta(beta[i,], beta[i-1,]))
    if (log(runif(1))>=acc.prob[i-1]){
      beta[i,-3] = beta[i-1,-3]
      if (i>burnin){
        browser()
        if(first){
          alfa = mlik.y1$alfa
          tau = mlik.y1$tau 
          first = F
        }else{
          alfa = alfa + mlik.y1$alfa
          tau = tau + mlik.y1$tau
        }
      }
      
    }else{
      if (i>burnin){
        browser()
        if(first){
          alfa = mlik.y2$alfa
          tau = mlik.y2$tau 
          first = F
        }else{
          alfa = alfa + mlik.y2$alfa
          tau = tau + mlik.y2$tau
        }
      }
      mlik.y1 = mlik.y2
    }
  }
  return(list(beta = beta, alfa = alfa/(N-burnin), tau = tau/(N-burnin), acc.prob))
}


set.seed(123)
df = sample.data.bi.linreg()
mod = bi.linreg.mcmc.w.inla(df)



