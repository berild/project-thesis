library(INLA)
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

fit.marginals <- function(ws,margs,len = 400){
  xmin <- min(margs[[1]])
  xmax <- max(margs[[1]])
  xx <- seq(xmin, xmax, len = len)
  tot_ws = sum(ws)
  marg = numeric(len)
  for (i in seq(nrow(margs[[1]]))){
    marg = marg + ws[i]/tot_ws*inla.dmarginal(xx, list(x = margs[[1]][i,], y = margs[[2]][i,]))
  }
  data.frame(x = xx, y = marg)
}

inla.w.is <- function(data, init, prior, d.prop, r.prop, fit.inla, n.prop,target.sd = 1){
  mlik = numeric(n.prop)
  eta = matrix(NA, ncol = length(init), nrow = n.prop)
  weight = numeric(n.prop)
  a.mean = init
  a.var = diag(5,length(init),length(init))
  pb <- txtProgressBar(min = 0, max = n.prop, style = 3)
  for (i in seq(n.prop)){
    setTxtProgressBar(pb, i)
    eta[i,] = r.prop(a.mean, sigma = a.var)
    mod = fit.inla(data,eta[i,])
    mlik[i] = mod$mlik
    weight[i] = mlik[i] + prior(eta[i,]) - d.prop(y = eta[i,], x = a.mean, sigma = a.var)
  }
  weight = exp(weight - max(weight))
  a.mean = c(sum(eta[,1]*weight)/sum(weight),sum(eta[,2]*weight)/sum(weight))
  margs = NA
  a.var = calc.var(a.mean,weight,eta)
  cat("Set mean =",a.mean,"; var =",a.var,"\nSampling...")
  eta = matrix(NA, ncol = length(init), nrow = n.prop)
  weight = numeric(n.prop)
  mlik = numeric(n.prop)
  for (i in seq(n.prop)){
    setTxtProgressBar(pb, i)
    eta[i,] = r.prop(a.mean, sigma = a.var)
    mod = fit.inla(data,eta[i,])
    margs = store.post(mod$dists,margs,i,n.prop)
    mlik[i] = mod$mlik
    weight[i] = mlik[i] + prior(eta[i,]) - d.prop(eta[i,], a.mean, sigma = a.var)
  }
  weight = exp(weight - max(weight))
  return(list(eta = eta,
              margs = margs,
              weight = weight))
}


sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

set.seed(1)
df = sample.linreg()
mod = inla.w.is(df,init = c(0,0), prior.beta, dq.beta, rq.beta, fit.inla, n.prop = 100)

fit.marginals <- function(ws,margs,len = 100){
  xmin <- min(margs[[1]])
  xmax <- max(margs[[1]])
  xx <- seq(xmin, xmax, len = len)
  marg = numeric(len)
  for (i in seq(nrow(margs[[1]]))){
    marg = marg + inla.dmarginal(xx, list(x = margs[[1]][i,], y = margs[[2]][i,]))
  }
  data.frame(x = xx, y = marg)
}


