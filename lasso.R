library(INLA)
library(tidyverse)
library(ISLR)
library(glmnet)
library(smoothmest)#Laplace distribution
library(mvtnorm)



fit.inla <- function(data, b) {
  data$oset <- data$x %*% matrix(b, ncol = 1)
  res <- inla(y ~ -1 + offset(oset), data = data)
  res <- inla.rerun(res)
  return(list(mlik = res$mlik[[1]], alfa = res$marginals.fixed[[1]], tau = res$marginals.hyperpar[[1]]))
}

prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))
  
  if(!log) { res <- exp(res) }
  
  return(res)
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
summary( fit.inla(d, b = rep(0, n.beta))$model )


x1 <-cbind(1,x)
ML.betas <- solve(t(x1)%*%x1)%*%t(x1)%*%y
ML.betas

stdev.samp <- .25 * solve(t(x)%*%x)


dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
}

lasso.mcmc.w.inla <- function(data, n.beta){
  N = 100000
  burnin = 500
  beta = matrix(data = NA,nrow = N, ncol = n.beta)
  beta[1,] = rep(0,n.beta)
  mod1 = fit.inla(data, beta[1,])
  tau = mod1$tau * 0
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    beta[i,] = rq.beta(beta[i-1,])
    mod2 = fit.inla(data,beta[i,])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.beta(beta[i,]) +
                   dq.beta(beta[i,],beta[i-1,]) - 
                   mod1$mlik -
                   prior.beta(beta[i-1,]) - 
                   dq.beta(beta[i-1,], beta[i,]))
    if (log(runif(1))>acc.prob[i-1]){
      beta[i,] = beta[i-1,]
      if (i >burnin){
        tau = tau + mod1$tau 
      }
    }else if (i > burnin){
      tau = tau + mod2$tau 
    }
  }
  return(list(beta = beta, 
              tau = tau/(N-burnin),
              acc.prob = min(exp(acc.prob),1)))
}

lasso.mcmc <- function(data, n.beta){
  N = 100000
  burnin = 500
  beta = matrix(data = NA,nrow = N, ncol = n.beta)
  beta[1,] = rep(0,n.beta)
  mod1 = fit.inla(data, beta[1,])
  tau = mod1$tau * 0
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    beta[i,] = rq.beta(beta[i-1,])
    mod2 = fit.inla(data,beta[i,])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.beta(beta[i,]) +
                   dq.beta(beta[i,],beta[i-1,]) - 
                   mod1$mlik -
                   prior.beta(beta[i-1,]) - 
                   dq.beta(beta[i-1,], beta[i,]))
    if (log(runif(1))>acc.prob[i-1]){
      beta[i,] = beta[i-1,]
      if (i >burnin){
        tau = tau + mod1$tau 
      }
    }else if (i > burnin){
      tau = tau + mod2$tau 
    }
  }
  return(list(beta = beta, 
              tau = tau/(N-burnin),
              acc.prob = sapply(exp(acc.prob),min,1)))
}

set.seed(123)
mod = lasso.mcmc.w.inla(d, n.beta)

save(mod, file = "lasso-mcmc-w-inla2.Rdata")
