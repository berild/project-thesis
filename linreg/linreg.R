library(INLA)
library(tidyverse)

rq.beta <- function(x, sigma = .75) {
  rnorm(length(x), mean = x, sd = sigma)
}

dq.beta <- function(y, x, sigma = .75, log =TRUE) {
  sum(dnorm(x, mean = y, sd = sigma, log = log))
}

prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}


fit.inla = function(data, beta){
  data$oset = beta[1]*data$x1 + beta[2]*data$x2
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              beta0 = res$marginals.fixed[[1]], 
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

dist.int <- function(x,n.step = 100){
   lx = min(x)
   ux = max(x)
   step = (ux- lx)/10
   seq(from = lx - step,to = ux + step, by = (ux - lx + 2*step)/(n.step-1))
}

moving.marginals <- function(new.marginal, avg.marginal, n){
  tmp = inla.dmarginal(avg.marginal$x, new.marginal, log = FALSE)
  (tmp + (n-1)*avg.marginal$y)/n
}

linreg.mcmc.w.inla <- function(data,n.samples = 100, n.burnin = 1, n.thin = 1){
  beta = matrix(data = NA,nrow = n.samples, ncol = 2)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  colnames(beta) = colnames(data[,-1])
  beta[1,] = c(0,0)
  mod.curr = fit.inla(data, beta = beta[1,])
  mlik[1] = mod.curr$mlik
  beta0 = data.frame(x = rep(0,100), y = rep(0,100))
  tau = data.frame(x = rep(0,100), y = rep(0,100))
  for (i in seq(2, n.samples)){
    beta.new = rq.beta(beta[i-1,])
    mod.new = fit.inla(data, beta = beta.new)
    lacc1 = mod.new$mlik + prior.beta(beta.new) + dq.beta(beta.new, beta[i-1,])
    lacc2 = mod.curr$mlik + prior.beta(beta[i-1,]) + dq.beta(beta[i-1,], beta.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      beta[i,] = beta.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
      
    }else{
      beta[i,] = beta[i-1,]
      mod.prev = mod.new
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i == n.burnin){
      browser()
      beta0$x = dist.int(x = mod.prev$beta0[,"x"])
      tau$x = dist.int(x = mod.prev$tau[,"x"])
    }else if(i > n.burnin){
      beta0$y = moving.marginals(new.marginal = mod.prev$beta0,
                                 avg.marginal = beta0,
                                 n = i-n)
      tau$y = moving.marginals(new.marginal = mod.prev$tau,
                                 avg.marginal = tau,
                                 n = i-n)
    }
  }
  return(list(beta = beta,
              acc.vec = acc.vec,
              mlik = mlik))
}


set.seed(1)
df = sample.linreg()
mod = linreg.mcmc.w.inla(df,n.samples = 1000,n.burnin = 100)
save(mod, file = "./linreg/linreg.Rdata")
mod_inla = inla(y~1 + x1 + x2,data = df)
save(mod_inla, file = "./linreg/linreg_INLA.Rdata")



