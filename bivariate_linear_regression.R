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
  return(list(mlik = res$mlik[1,1], model = res))
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
  beta = data.frame(beta_1 = 0, beta_2 = 0)
  mlik.y1 = fit.inla(df,beta[1,])
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in seq(2, N)){
    setTxtProgressBar(pb, i)
    beta[i,] = draw.prop.beta(beta[i-1,])
    mlik.y2 = fit.inla(df,beta[i,])
    acc.prob = mlik.y2$mlik + prior.beta(beta[i,]) +
      prob.prop.beta(beta[i-1,],beta[i,]) - mlik.y1$mlik -
      prior.beta(beta[i-1,]) - prob.prop.beta(beta[i,], beta[i-1,])
    if (log(runif(1))>=acc.prob){
      beta[i,] = beta[i-1,]
    }else{
      mlik.y1 = mlik.y2
    }
  }
  return(beta)
}


set.seed(123)
df = sample.data.bi.linreg()
beta = bi.linreg.mcmc.w.inla(df)
ggplot(beta,aes(x = seq(1,nrow(beta)))) + 
  geom_line(aes(y = beta_1, color = "beta_1"))+ 
  geom_line(aes(y = beta_2, color = "beta_2"))

beta2 = beta[seq(500,nrow(beta),10),]

ggplot(beta2,aes(x = beta_2)) + 
  geom_density()
