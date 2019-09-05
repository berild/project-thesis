library(INLA)
library(tidyverse)
library(ISLR)
library(glmnet)
library(smoothmest)#Laplace distribution
library(mvtnorm)


fit.inla <- function(data, b) {
  data$oset = data$x %*% b
  res = inla(y ~ -1 + offset(oset), data = data)
  res = inla.rerun(res)
  return(list(mlik = res$mlik[[1]],  
              tau = res$marginals.hyperpar[[1]]))
}

prior.beta <- function(x) {
  sum(log(ddoublex(x, mu = 0, lambda =  1/0.73)))
}


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

# proposal distribution of beta
dq.beta <- function(x, sigma = stdev.samp) {
  dmvnorm(x, mean = rep(0,length(x)), sigma = sigma, log = TRUE)
}

# sampling from the proposal of beta
rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = rep(0,length(x)), sigma = sigma))
}

lasso.mcmc.w.inla <- function(data, n.beta, stdev.samp){
  N = 10000
  burnin = 500
  beta = matrix(data = NA,nrow = N, ncol = n.beta)
  colnames(beta) = colnames(data$x)
  beta[1,] = rep(0,n.beta)
  mod1 = fit.inla(data, beta[1,])
  tau = mod1$tau * 0
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c(0)
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    beta[i,] = rq.beta(x = beta[i-1,],sigma = stdev.samp)
    mod2 = fit.inla(data,beta[i,])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.beta(beta[i,]) +
                   dq.beta(beta[i-1,], sigma = stdev.samp) - 
                   mod1$mlik -
                   prior.beta(beta[i-1,]) - 
                   dq.beta(beta[i,], sigma = stdev.samp))
    if (log(runif(1))>acc.prob[i]){
      beta[i,] = beta[i-1,]
    }else{ 
      mod1 = mod2
    }
    if (i > burnin){
      tau = tau + mod1$tau
    }
  }
  return(list(beta = beta, 
              tau = tau/(N-burnin),
              acc.prob = sapply(exp(acc.prob),min,1)))
}

set.seed(123)
mod = lasso.mcmc.w.inla(d, n.beta, stdev.samp)

save(mod, file = "./lasso/lasso-mcmc-w-inla-zero.Rdata")
