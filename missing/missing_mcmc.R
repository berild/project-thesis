library(mice)
library(INLA)

data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)


fit.inla <- function(data, x.mis) { 
  
  data$bmi[idx.mis] <- x.mis
  
  res <- inla(chl ~ 1 + bmi + age, data = data)
  
  return(list(mlik = res$mlik[[1]], 
              dists = list(beta0 = res$marginals.fixed[[1]], 
                           beta1 = res$marginals.fixed[[2]],
                           beta2 = res$marginals.fixed[[3]],
                           beta3 = res$marginals.fixed[[4]],
                           tau = res$marginals.hyperpar[[1]])))
}



dq.x.mis <- function(x, y, sigma = sqrt(10), log =TRUE) {
  res <- dnorm(y, mean = x, sd = sigma, log = log)
  
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

rq.x.mis <- function(x, sigma = sqrt(10) ) {
  rnorm(length(x), mean = x, sd = sigma)
}


prior.x.mis <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE), 
                        sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
  res = dnorm(x, mean = mu, sd= sigma, log = log)
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

r.beta <- function(data, beta, tau, sigma = 0.001){
  std = 
  uti = lapply(seq(length(beta)), function(x){sigma + tau*data$x})
}


r.tau <- function(data, beta, x.mis, idx.mis){
  data$bmi[idx.mis] <- x.mis
  shape = 1.5
  rate = 1/(5*10^(-5)) + 
    t(data$y - data$x%*%beta)%*%(data$y - data$x%*%beta)/2
  rgamma(1, shape = shape, rate = rate)
}

missing.mcmc <- function(data, n.mis,idx.mis, n.samples = 100, n.burnin = 5, n.thin = 1){
  chain = list(x.mis = matrix(NA,nrow = n.samples, ncol = n.mis),
               beta = matrix(NA,nrow = n.samples, ncol = 4),
               tau = numeric(n.samples),
               acc.vec = numeric(n.samples))
  mu = mean(data$bmi, na.rm = TRUE)
  sigma = 2*sd(data$bmi, na.rm = TRUE)
  chain$beta[1,] = rep(0,n.beta)
  chain$tau[1] = r.tau(data, chain$beta[1,])
  chain$acc.vec[1] = T
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
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