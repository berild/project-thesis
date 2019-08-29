library(mice)
library(INLA)

data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)


fit.inla <- function(data, x.mis) { 
  
  data$bmi[idx.mis] <- x.mis
  
  res <- inla(chl ~ 1 + bmi + age, data = data)
  
  return(list(mlik = res$mlik[1,1], model = res))
}



dq.beta <- function(x, y, sigma = sqrt(10), log =TRUE) {
  res <- dnorm(y, mean = x, sd = sigma, log = log)
  
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
  
  
}

rq.beta <- function(x, sigma = sqrt(10) ) {
  rnorm(length(x), mean = x, sd = sigma)
}


prior.beta <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE), 
                       sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
  res <- dnorm(x, mean = mu, sd= sigma, log = log)
  
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

inlamh.res <- INLAMH(d.mis, fit.inla, rep(mean(d.mis$bmi, na.rm = TRUE), n.mis),
                     rq.beta, dq.beta, prior.beta, 
                     n.sim = 10000, n.burnin = 500, n.thin = 10)

inlamh.res <- INLAMH(d, fit.inla, rep(0, n.beta), rq.beta, dq.beta, prior.beta,
                     n.sim = 10000, n.burnin = 500, n.thin = 10, verbose = TRUE)

missing.mcmc.w.inla <- function(data, n.mis){
  browser()
  N = 1000
  burnin = 500
  beta = matrix(data = NA,nrow = N, ncol = n.mis)
  beta[1,] = rep(data$bmi,n.mis)
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
  return(list(beta = beta, tau = tau/(N-burnin)))
}
