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

dist.init <- function(x,n.step = 100){
  lx = min(x)
  ux = max(x)
  step = (ux - lx)/(n.step-1)
  seq(from = lx - step, to = ux + step, by = (ux - lx)/(n.step-1))
}


moving.marginals <- function(new.marginal, avg.marginal, n){
  step = avg.marginal$x[2] - avg.marginal$x[1]
  new.x = avg.marginal$x
  new.y = avg.marginal$y
  if (max(new.marginal[,"x"])>max(new.x)){
    new.u = seq(from = max(new.x) + step,
                to = max(new.marginal[,"x"])+ step,
                by = step)
    new.x = c(new.x, new.u)
    new.y = c(new.y,rep(0,length(new.u)))
  }
  if (min(new.marginal[,"x"])<min(new.x)){
    new.l = seq(from = min(new.x) - step,
                to = min(new.marginal[,"x"]) - step,
                by = - step)
    new.x = c(rev(new.l), new.x)
    new.y = c(rep(0,length(new.l)),new.y)
  }
  avg.marginal = data.frame(x = new.x, y = new.y)
  tmp = inla.dmarginal(avg.marginal$x, new.marginal, log = FALSE)
  avg.marginal$y = (tmp + (n-1)*avg.marginal$y)/n
  avg.marginal
}

missing.mcmc.w.inla <- function(data, n.mis,idx.mis, n.samples = 100, n.burnin = 5, n.thin = 1){
  mu = mean(data$bmi, na.rm = TRUE)
  sigma = 2*sd(data$bmi, na.rm = TRUE)
  x.mis = matrix(data = NA,nrow = n.samples, ncol = n.mis)
  colnames(x.mis) = sprintf("Observation_%d",idx.mis)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  x.mis[1,] = rep(mu,n.mis)
  mod.curr = fit.inla(data,x.mis[1,])
  mlik[1] = mod.curr$mlik
  beta0 = data.frame(x = rep(0,102), y = rep(0,102))
  beta1 = data.frame(x = rep(0,102), y = rep(0,102))
  beta2 = data.frame(x = rep(0,102), y = rep(0,102))
  beta3 =data.frame(x = rep(0,102), y = rep(0,102))
  tau = data.frame(x = rep(0,102), y = rep(0,102))
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    x.mis.new = rq.x.mis(x.mis[i-1,])
    mod.new = fit.inla(data,x.mis[i,])
    lacc1 = mod.new$mlik + prior.x.mis(x.mis.new) + dq.x.mis(x.mis.new, x.mis[i-1,])
    lacc2 = mod.curr$mlik + prior.x.mis(x.mis[i-1,]) + dq.x.mis(x.mis[i-1,], x.mis.new)
    acc = acc = min(1,exp(lacc1 - lacc2))
    if (runif(1)>acc){
      x.mis[i,] = x.mis.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
    }else{ 
      x.mis[i,] = x.mis[i-1,]
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i == n.burnin){
      beta0$x = dist.init(x = mod.curr$beta0[,"x"])
      beta1$x = dist.init(x = mod.curr$beta1[,"x"])
      beta2$x = dist.init(x = mod.curr$beta2[,"x"])
      beta3$x = dist.init(x = mod.curr$beta3[,"x"])
      tau$x = dist.init(x = mod.curr$tau[,"x"])
    }else if(i > n.burnin){
      beta0 = moving.marginals(new.marginal = mod.curr$beta0,
                             avg.marginal = beta0,
                             n = i-n.burnin)
      beta1 = moving.marginals(new.marginal = mod.curr$beta1,
                             avg.marginal = beta1,
                             n = i-n.burnin)
      beta2 = moving.marginals(new.marginal = mod.curr$beta2,
                             avg.marginal = beta2,
                             n = i-n.burnin)
      beta3 = moving.marginals(new.marginal = mod.curr$beta3,
                             avg.marginal = beta3,
                             n = i-n.burnin)
      tau = moving.marginals(new.marginal = mod.curr$tau,
                             avg.marginal = tau,
                             n = i-n.burnin)
    }
    
  }
  return(list(x.mis = x.mis, 
              beta0 = beta0, 
              beta1 = beta1,
              beta2 = beta2,
              beta3 = beta3,
              tau = tau,
              mlik = mlik,
              acc.vec = acc.vec))
}

set.seed(123)
mod <- missing.mcmc.w.inla(d.mis, n.mis, idx.mis, n.samples = 100000, n.burnin = 500)
save(mod, file = "./missing/missing-mcmc-w-inla.Rdata")