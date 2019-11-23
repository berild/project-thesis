library(INLA)
library(tidyverse)
source("inla_w_mh.R")

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
  data$oset = data$x%*%beta
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]])))
}

sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

moving.marginals <- function(marg, post.marg, n){
  for (i in seq(length(post.marg))){
    tmp.post.marg = post.marg[[i]]
    tmp.marg = marg[[i]]
    step = tmp.post.marg[2,1] - tmp.post.marg[1,1]
    new.x = tmp.post.marg[,1]
    new.y = tmp.post.marg[,2]
    if (max(tmp.marg[,1])>max(new.x)){
      new.u = seq(from = max(new.x) + step,
                  to = max(tmp.marg[,1])+ step,
                  by = step)
      new.x = c(new.x, new.u)
      new.y = c(new.y,rep(0,length(new.u)))
    }
    if (min(tmp.marg[,1])<min(new.x)){
      new.l = seq(from = min(new.x) - step,
                  to = min(tmp.marg[,1]) - step,
                  by = - step)
      new.x = c(rev(new.l), new.x)
      new.y = c(rep(0,length(new.l)),new.y)
    }
    tmp.post.marg = data.frame(x = new.x, y = new.y)
    tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
    tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
    post.marg[[i]] = tmp.post.marg
  }
  post.marg
}

mcmc.w.inla <- function(data,n.samples = 100, n.burnin = 5, n.thin = 1){
  beta = matrix(data = NA,nrow = n.samples, ncol = 2)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  colnames(beta) = colnames(data[,-1])
  beta[1,] = c(0,0)
  mod.curr = fit.inla(data, beta = beta[1,])
  mlik[1] = mod.curr$mlik
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2, n.samples)){
    setTxtProgressBar(pb, i)
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
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i == n.burnin){
      post.marg = mod.curr$dists
    }else if(i > n.burnin){
      post.marg = moving.marginals(mod.curr$dists,
                                   post.marg,
                                   i-n.burnin+1)
    }
  }
  return(list(beta = beta,
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}


set.seed(1)
df = sample.linreg()
mod = inla.w.mcmc(df, c(0,0), prior.beta, rq.beta, dq.beta, fit.inla, n.samples = 100, n.burnin = 5, n.thin = 1)
save(mod, file = "./linreg/linreg.Rdata")
