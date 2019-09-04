library(mice)

data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)

dq.x.mis <- function(x, y, sigma = sqrt(0.001), log =TRUE) {
  sum(dnorm(y, mean = x, sd = sigma, log = log))
}

rq.x.mis <- function(x, sigma = sqrt(0.001)) {
  rnorm(length(x), mean = x, sd = sigma)
}


r.tau <- function(data, beta){
  shape = 1.5
  rate = 1/(5*10^(-5)) + t(data$y - data$x%*%beta)%*%(data$y - data$x%*%beta)
  rgamma(1, shape = shape, rate = rate)
}

d.y <- function(data,beta,tau,log = TRUE){
  dmvnorm(as.numeric(data$y),mean = as.numeric(data$x%*%beta), sigma = 1/tau*diag(length(data$y)),log = log)*prior.beta(beta)
}

prior.x.mis <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE), 
                        sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
  sum(dnorm(x, mean = mu, sd= sigma, log = log))
}

missing.mcmc <- function(data, n.mis){
  N = 10000
  burnin = 500
  mu = mean(data$bmi, na.rm = TRUE)
  sigma = 2*sd(data$bmi, na.rm = TRUE)
  x.mis = matrix(data = NA,nrow = N, ncol = n.mis)
  x.mis[1,] = rep(mean(data$bmi,na.rm = TRUE),n.mis)
  mod1 = fit.inla(data,x.mis[1,])
  alfa = mod1$alfa * 0
  beta = mod1$beta * 0 
  tau = mod1$tau * 0
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c()
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    x.mis[i,] = rq.x.mis(x.mis[i-1,])
    mod2 = fit.inla(data,x.mis[i,])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.x.mis(x.mis[i,],mu = mu, sigma = sigma) +
                   dq.x.mis(x.mis[i,],x.mis[i-1,]) - 
                   mod1$mlik -
                   prior.x.mis(x.mis[i-1,], mu = mu, sigma = sigma) - 
                   dq.x.mis(x.mis[i-1,], x.mis[i,]))
    if (log(runif(1))>acc.prob[i-1]){
      x.mis[i,] = x.mis[i-1,]
      if (i >burnin){
        alfa = alfa +  mod1$alfa
        beta = beta + mod1$beta
        tau = tau + mod1$tau 
      }
    }else if (i > burnin){
      alfa = alfa +  mod2$alfa
      beta = beta + mod2$beta
      tau = tau + mod2$tau 
    }
  }
  return(list(x.mis = x.mis, 
              alfa = alfa/(N-burnin), 
              beta = beta/(N-burnin), 
              tau = tau/(N-burnin),
              acc.prob = sapply(exp(acc.prob),min,1)))
}

mod <- missing.mcmc(d.mis, n.mis)

save(mod, file = "missing-mcmc.Rdata")