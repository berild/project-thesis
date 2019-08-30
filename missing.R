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
              alfa = res$marginals.fixed[[1]], 
              beta = res$marginals.fixed[[2]], 
              tau = res$marginals.hyperpar[[1]]))
}



dq.x.mis <- function(x, y, sigma = sqrt(10), log =TRUE) {
  sum(dnorm(y, mean = x, sd = sigma, log = log))
}

rq.x.mis <- function(x, sigma = sqrt(10) ) {
  rnorm(length(x), mean = x, sd = sigma)
}


prior.x.mis <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE), 
                       sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
  sum(dnorm(x, mean = mu, sd= sigma, log = log))
}

missing.mcmc.w.inla <- function(data, n.mis){
  N = 10000
  burnin = 500
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
                   prior.x.mis(x.mis[i,]) +
                   dq.x.mis(x.mis[i,],x.mis[i-1,]) - 
                   mod1$mlik -
                   prior.x.mis(x.mis[i-1,]) - 
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

mod <- missing.mcmc.w.inla(d.mis, n.mis)

save(mod, file = "missing-mcmc-w-inla.Rdata")