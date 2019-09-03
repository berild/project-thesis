library(INLA)
library(tidyverse)

draw.prop.beta <- function(b,sigma = 1/0.75){
  rnorm(ncol(b),mean = as.numeric(b), sd = sigma)
}

prob.prop.beta <- function(b1,b2, sigma = 1/0.75, log = TRUE){
  sum(dnorm(as.numeric(b1),mean = as.numeric(b2), sd = sigma, log = log))
}

prior.beta <- function(b, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(as.numeric(b), mean = 0, sd= sigma, log = log))
}


fit.inla <- function(data,b){
  data$oset = b[1] * data$x1 +b[2]*data$x2
  res = inla(y ~1+offset(oset), data = data)
  return(list(mlik = res$mlik[[1]], alfa = res$marginals.fixed[[1]], tau = res$marginals.hyperpar[[1]], mod = res))
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
  N = 1000
  burnin = 500
  beta = data.frame(beta_1 = rep(NA,N),beta_2 = rep(NA,N))
  beta[1,] = c(0,0)
  mlik.y1 = fit.inla(df,as.numeric(beta[1,]))
  first = T
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  browser()
  acc.prob = c()
  for (i in seq(2, N)){
    setTxtProgressBar(pb, i)
    beta[i,] = draw.prop.beta(beta[i-1,])
    mlik.y2 = fit.inla(df,as.numeric(beta[i,]))
    acc.prob = c(acc.prob,
                 mlik.y2$mlik + 
                   prior.beta(beta[i,]) +
                   prob.prop.beta(beta[i-1,],beta[i,]) - 
                   mlik.y1$mlik -
                   prior.beta(beta[i-1,]) - 
                   prob.prop.beta(beta[i,], beta[i-1,]))
    if (log(runif(1))>=acc.prob[i-1]){
      beta[i,] = beta[i-1,]
      if (i>burnin){
        if(first){
          alfa = mlik.y1$alfa
          tau = mlik.y1$tau 
          first = F
        }else{
          alfa = alfa + mlik.y1$alfa
          tau = tau + mlik.y1$tau
        }
      }
      
    }else{
      if (i>burnin){
        if(first){
          alfa = mlik.y2$alfa
          tau = mlik.y2$tau 
          first = F
        }else{
          alfa = alfa + mlik.y2$alfa
          tau = tau + mlik.y2$tau
        }
      }
      mlik.y1 = mlik.y2
    }
  }
  beta = cbind(beta, is_burnin = c(rep(T,burnin)),rep(F,N-burnin))
  return(list(beta = beta, alfa = alfa/(N-burnin), tau = tau/(N-burnin), acc.prob = acc.prob))
}


set.seed(123)
df = sample.data.bi.linreg()
mod = bi.linreg.mcmc.w.inla(df)
mod.inla = inla(y~1 + x1 + x2, data = df)



ggplot(as.data.frame(mod$alfa), aes(x = x, y = y)) + 
  geom_line()

ggplot(as.data.frame(mod$tau), aes(x = x, y = y)) + 
  geom_line()

ggplot(as.data.frame(mod$beta)) + 
  geom_density(aes(beta_1, fill = "beta_1"),alpha = 0.5) + 
  geom_density(aes(beta_2, fill = "beta_2"),alpha = 0.5) + 
  labs(fill = "") 

ggplot(as.data.frame(mod$beta),aes(x = seq(1,length(beta_1)))) + 
  geom_line(aes(y = beta_1, color = "beta_1")) + 
  geom_line(aes(y = beta_2, color = "beta_2")) + 
  labs(color = "")
