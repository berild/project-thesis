library(tidyverse)
library(ISLR)
library(glmnet)
library(smoothmest) #Laplace distribution
library(mvtnorm)



prior.beta <- function(x, mu = 0, lambda = 1/0.73, log = TRUE) {
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

r.tau <- function(data, beta){
  shape = 1.5
  rate = 1/(5*10^(-5)) + t(data$y - data$x%*%beta)%*%(data$y - data$x%*%beta)
  rgamma(1, shape = shape, rate = rate)
}

d.y <- function(data,beta,tau,log = TRUE){
  dmvnorm(as.numeric(data$y),mean = as.numeric(data$x%*%beta), sigma = 1/tau*diag(length(data$y)),log = log)*prior.beta(beta)
}

lasso.mcmc <- function(data, n.beta){
  N = 100000
  beta = matrix(data = NA,nrow = N, ncol = n.beta)
  beta[1,] = rep(0,n.beta)
  tau = c()
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    tau = c(tau, r.tau(data,beta[i-1,]))
    beta[i,] = rq.beta(beta[i-1,])
    acc.prob = c(acc.prob,
                 d.y(data,beta[i,],tau[i-1]) + 
                   dq.beta(beta[i,],beta[i-1,]) - 
                   d.y(data,beta[i-1,],tau[i-1]) - 
                   dq.beta(beta[i-1,], beta[i,]))
    if (log(runif(1))>acc.prob[i-1]){
      beta[i,] = beta[i-1,]
    }
  }
  return(list(beta = beta, 
              tau = tau,
              acc.prob = sapply(exp(acc.prob),min,1)))
}

set.seed(123)
  mod = lasso.mcmc(d, n.beta)

save(mod, file = "lasso-mcmc.Rdata")
