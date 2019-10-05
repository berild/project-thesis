library(tidyverse)
library(ISLR)
library(glmnet)
library(smoothmest) #Laplace distribution
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


r.tau <- function(data, beta){
  shape = 1.5
  rate = 5*10^(-5) + t(data$y - data$x%*%beta)%*%(data$y - data$x%*%beta)/2
  rgamma(1, shape = shape, rate = rate)
}

d.beta <- function(data,beta,tau,log = TRUE){
  dmvnorm(as.numeric(data$y),
          mean = as.numeric(data$x%*%beta), 
          sigma = 1/tau*diag(length(data$y)),
          log = log) + prior.beta(beta)
}

lasso.mcmc <- function(data, n.beta, stdev.samp, n.samples = 100, n.burnin = 5, n.thin = 1){
  chain = list(beta = matrix(NA,nrow = n.samples, ncol = n.beta),
               tau = numeric(n.samples),
               acc.vec = numeric(n.samples))
  chain$beta[1,] = rep(0,n.beta)
  chain$tau[1] = r.tau(data, chain$beta[1,])
  chain$acc.vec[1] = T
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    beta.new = rq.beta(chain$beta[i-1,],stdev.samp)
    tau.new = r.tau(data,beta.new)
    lacc1 = d.beta(data,beta.new,tau.new) + dq.beta(beta.new,chain$beta[i-1,],stdev.samp)
    lacc2 = d.beta(data,chain$beta[i-1,],tau.new) + dq.beta(chain$beta[i-1,],beta.new,stdev.samp)
    if (runif(1)<exp(lacc1 - lacc2)){
      chain$beta[i,] = beta.new
      chain$tau[i] = tau.new
      chain$acc.vec[i] = T
    }else{
      chain$beta[i,] = chain$beta[i-1,]
      chain$tau[i] = chain$tau[i-1]
      chain$acc.vec[i] = F
    }
  }
  return(chain)
}


set.seed(123)
mod = lasso.mcmc(d, n.beta,stdev.samp,n.samples = 1000)



mean(mod$beta[,1])
save(mod, file = "./lasso/lasso-mcmc.Rdata")
