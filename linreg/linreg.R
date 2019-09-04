library(INLA)
library(tidyverse)

draw.prop.beta <- function(b,sigma = 1/0.75){
  rnorm(length(b),mean = b, sd = sigma)
}

prob.prop.beta <- function(b1,b2, sigma = 1/0.75, log = TRUE){
  sum(dnorm(b1,mean = b2, sd = sigma, log = log))
}

prior.beta <- function(b, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(b, mean = 0, sd= sigma, log = log))
}


fit.inla <- function(data,b){
  data$oset = b[1] * data$x1 + b[2]*data$x2
  res = inla(y ~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1,1], 
              alfa = res$marginals.fixed[[1]], 
              tau = res$marginals.hyperpar[[1]]))
}


sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2
  df = data.frame(y = y, x1 = x1, x2 = x2)
}


linreg.mcmc.w.inla <- function(data){
  N = 10000
  burnin = 500
  beta = matrix(data = NA,nrow = N, ncol = 2)
  colnames(beta) = colnames(data[,-1])
  beta[1,] = c(0,0)
  mod1 = fit.inla(data,beta[1,])
  alfa = mod1$alfa*0
  tau = mod1$tau * 0
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2, N)){
    setTxtProgressBar(pb, i)
    beta[i,] = draw.prop.beta(beta[i-1,])
    mod2 = fit.inla(data,beta[i,])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.beta(beta[i,]) +
                   prob.prop.beta(beta[i-1,],beta[i,]) - 
                   mod1$mlik -
                   prior.beta(beta[i-1,]) - 
                   prob.prop.beta(beta[i,], beta[i-1,]))
    if (log(runif(1))>=acc.prob[i-1]){
      beta[i,] = beta[i-1,]
      
    }else{
      mod1 = mod2
    }
    if (i > burnin){
      alfa = alfa + mod1$alfa
      tau = tau + mod1$tau
    }
  }
  return(list(alfa = alfa/(N-burnin), 
              beta = beta,
              tau = tau/(N-burnin), 
              acc.prob = sapply(exp(acc.prob),min,1)))
}


set.seed(123)
df = sample.linreg()
mod = linreg.mcmc.w.inla(df)
save(mod, file = "./linreg/linreg.Rdata")
mod_inla = inla(y~1 + x1 + x2,data = df)
save(mod_inla, file = "./linreg/linreg_INLA.Rdata")