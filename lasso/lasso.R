library(INLA)
library(tidyverse)
library(ISLR)
library(glmnet)
library(smoothmest)
library(mvtnorm)


data(Hitters)

Hitters <- na.omit(Hitters)

x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5]
y <- Hitters$Salary


y <- scale(y)
x <- scale(x)
d <- list(y = y, x = x)
n.beta <- ncol(d$x)

# finding inverse of the precision
stdev.samp <- .25 * solve(t(x)%*%x)


fit.inla <- function(data, b) {
  data$oset = data$x %*% b
  res = inla(y ~ -1 + offset(oset), data = data)
  res = inla.rerun(res)
  return(list(mlik = res$mlik[[1]],  
              tau = res$marginals.hyperpar[[1]]))
}

prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))
  
  if(!log) { res <- exp(res) }
  
  return(res)
}

dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
}

dist.init <- function(x,n.step = 100){
  lx = min(x)
  ux = max(x)
  step = (ux - lx)/(n.step-1)
  seq(from = lx - step, to = ux + step, by = (ux - lx)/(n.step-1))
}


moving.marginals <- function(new.marginal, avg.marginal, n){
  step = avg.marginal$x[2] - avg.marginal$x[1]
  new.x = avg.marginal$x
  new.y = avg.marginal$y
  if (max(new.marginal[,"x"])>max(new.x)){
    new.u = seq(from = max(new.x) + step,
                to = max(new.marginal[,"x"])+ step,
                by = step)
    new.x = c(new.x, new.u)
    new.y = c(new.y,rep(0,length(new.u)))
  }
  if (min(new.marginal[,"x"])<min(new.x)){
    new.l = seq(from = min(new.x) - step,
                to = min(new.marginal[,"x"]) - step,
                by = - step)
    new.x = c(rev(new.l), new.x)
    new.y = c(rep(0,length(new.l)),new.y)
  }
  avg.marginal = data.frame(x = new.x, y = new.y)
  tmp = inla.dmarginal(avg.marginal$x, new.marginal, log = FALSE)
  avg.marginal$y = (tmp + (n-1)*avg.marginal$y)/n
  avg.marginal
}


lasso.mcmc.w.inla <- function(data, n.beta, stdev.samp, n.samples = 100, n.burnin = 5, n.thin = 1){
  beta = matrix(data = NA,nrow = n.samples, ncol = n.beta)
  colnames(beta) = colnames(data$x)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  beta[1,] = rep(0,n.beta)
  mod.curr = fit.inla(data, beta[1,])
  mlik[1] = mod.curr$mlik
  tau = data.frame(x = rep(0,102), y = rep(0,102))
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    beta.new = rq.beta(beta[i-1,],sigma = stdev.samp)
    mod.new = fit.inla(data,beta.new)
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
      tau$x = dist.init(x = mod.curr$tau[,"x"])
    }else if(i > n.burnin){
      tau = moving.marginals(new.marginal = mod.curr$tau,
                             avg.marginal = tau,
                             n = i-n.burnin)
    }
  }
  return(list(beta = beta, 
              tau = tau,
              acc.vec = acc.vec))
}

set.seed(123)
mod = lasso.mcmc.w.inla(d, n.beta, stdev.samp, n.samples = 100000, n.burnin = 500)
save(mod, file = "./lasso/lasso-mcmc-w-inla.Rdata")
