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

moving.marginals <- function(marg, post.marg, n){
  for (i in seq(length(post.marg))){
    tmp.post.marg = post.marg[[i]]
    tmp.marg = marg[[i]]
    step = tmp.post.marg[2,1] - tmp.post.marg[1,1]
    new.x = tmp.post.marg[,1]
    new.y = tmp.post.marg[,2]
    if (max(tmp.marg[,1])>max(new.x)){
      new.u = seq(from = max(new.x) + step,
                  to = max(tmp.marg[,1])+ step,
                  by = step)
      new.x = c(new.x, new.u)
      new.y = c(new.y,rep(0,length(new.u)))
    }
    if (min(tmp.marg[,1])<min(new.x)){
      new.l = seq(from = min(new.x) - step,
                  to = min(tmp.marg[,1]) - step,
                  by = - step)
      new.x = c(rev(new.l), new.x)
      new.y = c(rep(0,length(new.l)),new.y)
    }
    tmp.post.marg = data.frame(x = new.x, y = new.y)
    tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
    tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
    post.marg[[i]] = tmp.post.marg
  }
  post.marg
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
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    x.mis.new = rq.x.mis(x.mis[i-1,])
    mod.new = fit.inla(data,x.mis[i,])
    lacc1 = mod.new$mlik + prior.x.mis(x.mis.new) + dq.x.mis(x.mis.new, x.mis[i-1,])
    lacc2 = mod.curr$mlik + prior.x.mis(x.mis[i-1,]) + dq.x.mis(x.mis[i-1,], x.mis.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1)<acc){
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
      post.marg = mod.curr$dists
    }else if(i > n.burnin){
      post.marg = moving.marginals(mod.curr$dists,
                                   post.marg,
                                   i-n.burnin+1)
    }
    
  }
  return(list(x.mis = x.mis, 
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}

set.seed(123)
mod <- missing.mcmc.w.inla(d.mis, n.mis, idx.mis, n.samples = 100000, n.burnin = 500)
save(mod, file = "./missing/missing-mcmc-w-inla.Rdata")