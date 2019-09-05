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
              beta0 = res$marginals.fixed[[1]], 
              beta1 = res$marginals.fixed[[2]],
              beta2 = res$marginals.fixed[[3]],
              beta3 = res$marginals.fixed[[4]],
              tau = res$marginals.hyperpar[[1]]))
}



dq.x.mis <- function(x, y) {
  sum(dnorm(y, mean = x, sd = sqrt(10), log = TRUE))
}

rq.x.mis <- function(x) {
  rnorm(length(x), mean = x, sd = sqrt(10))
}


prior.x.mis <- function(x, mu, sigma) {
  sum(dnorm(x, mean = mu, sd= sigma, log = TRUE))
}

missing.mcmc.w.inla <- function(data, n.mis,idx.mis){
  N = 100000
  burnin = 500
  mu = mean(data$bmi, na.rm = TRUE)
  sigma = 2*sd(data$bmi, na.rm = TRUE)
  x.mis = matrix(data = NA,nrow = N, ncol = n.mis)
  colnames(x.mis) = sprintf("Observation_%d",idx.mis)
  x.mis[1,] = rep(mu,n.mis)
  mod1 = fit.inla(data,x.mis[1,])
  beta0 = mod1$alfa * 0
  beta1 = mod1$beta * 0 
  beta2 = mod1$beta * 0 
  beta3 = mod1$beta * 0 
  tau = mod1$tau * 0
  pb = txtProgressBar(min = 0, max = N, style = 3)
  acc.prob = c(0)
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
    if (log(runif(1))>acc.prob[i]){
      x.mis[i,] = x.mis[i-1,]
    }else{ 
      mod1 = mod2
    }
    if (i >burnin){
      beta0 = beta0 +  mod1$beta0
      beta1 = beta1 + mod1$beta1
      beta2 = beta2 + mod1$beta2
      beta3 = beta3 + mod1$beta3
      tau = tau + mod1$tau 
    }
  }
  return(list(x.mis = x.mis, 
              beta0 = beta0/(N-burnin), 
              beta1 = beta1/(N-burnin),
              beta2 = beta2/(N-burnin),
              beta3 = beta3/(N-burnin),
              tau = tau/(N-burnin),
              acc.prob = sapply(exp(acc.prob),min,1)))
}

set.seed(123)
mod <- missing.mcmc.w.inla(d.mis, n.mis, idx.mis)

save(mod, file = "./missing/missing-mcmc-w-inla-large.Rdata")