library(INLA)
library(tidyverse)

rq.beta <- function(x=c(0,0), sigma = 5) {
  rnorm(length(x), mean = x, sd = sigma)
}

dq.beta <- function(y, x, sigma = 5, log =TRUE) {
  sum(dnorm(y, mean = x, sd = sigma, log = log))
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


inla.w.is <- function(data, init, prior, d.prop, r.prop, n.prop,target.sd = 1){
  mlik = numeric(n.prop)
  eta = matrix(NA, ncol = length(init), nrow = n.prop)
  weight = numeric(n.prop)
  post.marg = NA
  a.mean = init
  a.sd = 5
  m.weight = 0
  n = 3
  pb <- txtProgressBar(min = 0, max = n.prop, style = 3)
  for (i in seq(n.prop)){
    setTxtProgressBar(pb, i)
    eta[i,] = r.prop(a.mean, sigma = a.sd)
    mod = fit.inla(data,eta[i,])
    mlik[i] = mod$mlik
    weight[i] = exp(mlik[i]) #+ prior(eta[i,]) - d.prop(eta[i,], a.mean, sigma = a.sd))
    #a.sd = adaptive.sd2(a.sd,a.mean,weight[i], m.weight, eta[i,], n)
    #a.sd = adaptive.sd(a.sd, weight[i], m.weight,target.sd)
    a.sd = adaptive.sd3(a.sd, a.mean, weight[i], m.weight, eta[i,])
    a.mean = adaptive.mean(a.mean, weight[i], m.weight, eta[i,])
    post.marg = self.norm.imp.samp(mod$dists,post.marg,weight[i],m.weight)
    m.weight = m.weight + weight[i]
    n = n+1
  }
  return(list(eta = eta,
              post.marg = post.marg,
              weight = weight,
              a.mean = a.mean,
              a.sd = a.sd))
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
mod = inla.w.is(df,init = c(0,0), prior.beta, dq.beta, rq.beta, n.prop = 500)

ggplot(as.data.frame(mod$post.marg$intercept)) + 
  geom_line(aes(x = x, y = y))


ggplot(as.data.frame(mod$post.marg$tau)) + 
  geom_line(aes(x = x, y = y))

ggplot(as.data.frame(mod$eta)) + 
  geom_histogram(aes(x = V1),bins = 30)

ggplot(as.data.frame(mod$eta)) + 
  geom_histogram(aes(x = V2),bins = 30)

mod$eta
mod$weight
mod$a.mean
mod$a.sd

