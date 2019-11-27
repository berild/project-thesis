library(INLA)
library(tidyverse)
library(mvtnorm)

rq.beta <- function(x=c(0,0), sigma = diag(5,2,2)) {
  rmvnorm(1,mean=x,sigma = sigma)
}

dq.beta <- function(y, x, sigma = diag(5,2,2), log =TRUE) {
  dmvnorm(y,mean = x, sigma = sigma,log = log)
}

prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}


fit.inla = function(data, beta){
  data$oset = data$x%*%beta
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]])))
}


self.norm.imp.samp <- function(marg, post.marg = NA, weight, m.weight){
  for (i in seq(length(marg))){
    if(length(post.marg)!=length(marg)){
      return(marg)
    }else{
      tmp.post.marg = post.marg[[i]]
      tmp.marg = marg[[i]]
      min.x = min(c(tmp.post.marg[1,1],tmp.marg[1,1]))
      max.x = max(c(tmp.post.marg[length(tmp.post.marg[,1]),1], tmp.marg[length(tmp.marg[,1]),1]))
      new.x = seq(min.x,max.x,length.out = 100)
      new.y = inla.dmarginal(new.x, tmp.post.marg, log = FALSE)
      # step1 = tmp.post.marg[2,1] - tmp.post.marg[1,1]
      # step2 = tmp.marg[2,1] - tmp.marg[1,1]
      # if (step2 < step1){
      #   new.x = seq(tmp.post.marg[1,1],tmp.post.marg[length(tmp.post.marg[,1]),1],step2)
      #   new.y = inla.dmarginal(new.x,tmp.post.marg, log = FALSE)
      #   tmp.post.marg = data.frame(x = new.x,y = new.y)
      #   step = step2
      # }else{
      #   step = step1
      # }
      # new.x = tmp.post.marg[,1]
      # new.y = tmp.post.marg[,2]
      # if (max(tmp.marg[,1])>max(new.x)){
      #   new.u = seq(from = max(new.x) + step,
      #               to = max(tmp.marg[,1])+ step,
      #               by = step)
      #   new.x = c(new.x, new.u)
      #   new.y = c(new.y,rep(0,length(new.u)))
      # }
      # if (min(tmp.marg[,1])<min(new.x)){
      #   new.l = seq(from = min(new.x) - step,
      #               to = min(tmp.marg[,1]) - step,
      #               by = - step)
      #   new.x = c(rev(new.l), new.x)
      #   new.y = c(rep(0,length(new.l)),new.y)
      # }
      tmp.post.marg = data.frame(x = new.x, y = new.y)
      tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
      tmp.post.marg[,2] = tmp*weight/(m.weight + weight) + tmp.post.marg[,2]*m.weight/(m.weight + weight)
      post.marg[[i]] = tmp.post.marg
    }
    return(post.marg)
  }
}

adaptive.mean <- function(a.mean, new.weight, m.weight, new.eta){
  a.mean*m.weight/(m.weight + new.weight) + new.weight/(m.weight + new.weight)*new.eta 
}

adaptive.sd <- function(a.sd, new.weight, m.weight, target.sd){
  a.sd*m.weight/(m.weight + new.weight) + new.weight/(m.weight + new.weight)*target.sd
}

adaptive.sd2 <- function(a.sd, a.mean, new.weight, m.weight, new.eta, n){
  sqrt((n-2)/(n-1)*a.sd^2 + 1/n*t(new.weight/(new.weight + m.weight)*new.eta - m.weight/(new.weight + m.weight)*a.mean)%*%
         (new.weight/(new.weight + m.weight)*new.eta - m.weight/(new.weight + m.weight)*a.mean))
}

adaptive.sd3 <-function(a.sd, a.mean, new.weight, m.weight, new.eta){
  sqrt(m.weight^2/(new.weight + m.weight)^2*a.sd^2 + new.weight^2/(new.weight + m.weight)^2*t(new.eta - a.mean)%*%(new.eta - a.mean))
}

calc.var <- function(a.mean,weight,eta){
  new.var = matrix(NA,ncol = length(a.mean), nrow = length(a.mean))
  for (i in seq(length(a.mean))){
    for (j in seq(i,length(a.mean))){
      new.var[i,j] = new.var[j,i] = sum(weight*weight*(eta[,i]-a.mean[i])*(eta[,j]-a.mean[j]))/(sum(weight))^2
    }
  }
  return(new.var)
}

store.post <- function(marg,margs,j,n.prop){
  if (length(marg)!=length(margs)){
    margs = marg
    for (i in seq(length(marg))){
      margs[[i]] = list(x = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,1])),
                        y = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,2])))
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }else{
    for (i in seq(length(marg))){
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }
}



is.w.inla <- function(data, init, prior, d.prop, r.prop, N_0 = 200, N = 400){
  mlik = numeric(N_0)
  eta = matrix(NA, ncol = length(init$mu), nrow = N_0)
  weight = numeric(N_0)
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), 2))) 
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  pb <- txtProgressBar(min = 0, max = N+N_0, style = 3)
  for (i in seq(N_0)){
    setTxtProgressBar(pb, i)
    eta[i,] = r.prop(theta$a.mu[1,], sigma = theta$a.cov[,,1])
    mod = fit.inla(data,eta[i,])
    mlik[i] = mod$mlik
    weight[i] = mlik[i] + prior(eta[i,]) - d.prop(y = eta[i,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1])
  }
  weight = exp(weight - max(weight))
  theta = calc.theta(theta,weight,eta,nrow(eta),2)
  margs = NA
  eta = matrix(NA, ncol = length(init$mu), nrow = N)
  weight = numeric(N)
  mlik = numeric(N)
  for (i in seq(N)){
    setTxtProgressBar(pb, i+N_0)
    eta[i,] = r.prop(theta$a.mu[2,], sigma = theta$a.cov[,,2])
    mod = fit.inla(data,eta[i,])
    margs = store.post(mod$dists,margs,i,N)
    mlik[i] = mod$mlik
    weight[i] = mlik[i] + prior(eta[i,]) - d.prop(y = eta[i,], x = theta$a.mu[2,], sigma = theta$a.cov[,,2])
  }
  weight = exp(weight - max(weight))
  eta_kern = kde2d.weighted(x = eta[,1], y = eta[,2], w = weight/(sum(weight)), n = 100, lims = c(1,3,-3,-1))
  eta_kern = data.frame(expand.grid(x=eta_kern$x, y=eta_kern$y), z=as.vector(eta_kern$z))
  return(list(eta = eta,
              theta = theta,
              eta_kern = eta_kern,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight))
}
