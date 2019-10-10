library(tidyverse)
library(mvtnorm)

rq.beta <- function(x, sigma = .75) {
  rnorm(length(x), mean = x, sd = sigma)
}

dq.beta <- function(y, x, sigma = .75, log =TRUE) {
  sum(dnorm(x, mean = y, sd = sigma, log = log))
}

r.tau <- function(y,x,beta,alpha,a = 1,b = 5e-05){
  rgamma(1, shape = a + 0.5, rate = b + t(y- alpha - x%*%beta)%*%(y- alpha - x%*%beta)/2)
}

r.alpha <- function(y,x,beta,tau){
  rnorm(1, mean = sum(y - x%*%beta)/length(y), sd = 1 / sqrt(length(y)*tau))
}


prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}

d.beta <- function(y,x,alpha,beta,tau,log = TRUE){
  return(dmvnorm(y, mean = alpha + x%*%beta, sigma = 1/tau*diag(length(y)), log = log) +
           prior.beta(beta))
}

sample.linreg <- function(){
  n = 100
  x = matrix(runif(2*n),nrow = n)
  list(y = 3 + 2*x[,1] -2*x[,2] + rnorm(n), x = x )
}


# implement gibbs sampling here
# mcmc doens't work well with this model
# compare distributions
# deliver results to sara
# start work on spatial econometric models
# setup text
# then maybe add more

linreg.mcmc <- function(data,n.samples=100,n.burnin=5,n.thin = 1){
  chain = list(alpha = numeric(n.samples),
               beta = matrix(NA,nrow = n.samples, ncol = 2),
               tau = numeric(n.samples),
               acc.vec = numeric(n.samples))
  chain$beta[1,] = c(0,0)
  chain$alpha[1] = mean(data$y)
  chain$tau[1] = r.tau(data$y,data$x,chain$beta[1,],chain$alpha[1])
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    beta.new = rq.beta(chain$beta[i-1,])
    tau.new = r.tau(data$y,data$x,beta.new,chain$alpha[i-1])
    alpha.new = r.alpha(data$y,data$x,beta.new,tau.new)
    lacc1 = d.beta(data$y,data$x,alpha.new,beta.new,tau.new) + dq.beta(beta.new,chain$beta[i-1,])
    lacc2 = d.beta(data$y,data$x, chain$alpha[i-1],chain$beta[i-1,],chain$tau[i-1]) + dq.beta(chain$beta[i-1,],beta.new)
    if (runif(1)<exp(lacc1 - lacc2)){
      chain$alpha[i] = alpha.new
      chain$beta[i,] = beta.new
      chain$tau[i] = tau.new
      chain$acc.vec[i] = T
    }else{
      chain$alpha[i] = chain$alpha[i-1]
      chain$beta[i, ] = chain$beta[i-1,]
      chain$tau[i] = chain$tau[i-1]
      chain$acc.vec[i] = F
    }
    
  }
  return(chain)
}


set.seed(1)
samps = sample.linreg()
mod = linreg.mcmc(samps,n.samples = 10000)

res = data.frame(step = seq(10000),alpha = mod$alpha,beta1 = mod$beta[,1],beta2 = mod$beta[,2],tau = mod$tau)
ggplot(res, aes(x = step)) + 
  geom_line(aes(y = alpha))


