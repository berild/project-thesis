library(mvtnorm)
library(INLA)
library(MASS)

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

calc.cov <- function(a.mean,weight,eta){
  new.var = matrix(NA,ncol = length(a.mean), nrow = length(a.mean))
  for (i in seq(length(a.mean))){
    for (j in seq(i,length(a.mean))){
      new.var[i,j] = new.var[j,i] = sum(weight*(eta[,i]-a.mean[i])*(eta[,j]-a.mean[j]))/(sum(weight))
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

calc.delta <- function(N_t,eta,theta,t,d.prop){
  tmp = 0
  for (l in seq(t)){
    tmp = tmp + N_t[l]*d.prop(y = eta, x = theta$a.mu[l+1,], sigma = theta$a.cov[,,l+1], log = FALSE)
  }
  return(tmp)
}

update.delta.weight <- function(delta,weight,N_t,eta,theta,t,mlik,prior,d.prop){
  i_tmp = 0
  N_tmp = sum(N_t[1:(t+1)])
  for (l in seq(t)){
    for (i in seq(N_t[l])){
      i_tmp = i_tmp + 1
      delta[i_tmp] = delta[i_tmp] + N_t[l]*d.prop(y = eta[i_tmp,], x = theta$a.mu[t+1,], sigma = theta$a.cov[,,t+1], log = FALSE)
      weight[i_tmp] = exp(mlik[i_tmp] + prior(eta[i_tmp,]))/(delta[i_tmp]/N_tmp)
    }
  }
  return(list(
    delta = delta,
    weight = weight
  ))
}

inla.w.is <- function(data, init, prior, d.prop, r.prop, N_t = rep(20,20)){
  N_0 = round(sum(N_t)/2)
  N_tot = N_0 + sum(N_t)
  mlik = numeric(N_tot)
  eta = matrix(NA, ncol = length(init), nrow = N_tot)
  delta = numeric(N_tot)
  weight = numeric(N_tot)
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = length(N_t) + 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), length(N_t) +2))) 
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  
  # initialization process 
  i_tot = 0
  pb <- txtProgressBar(min = 0, max = N_tot, style = 3)
  margs = NA
  for (i in seq(N_0)){
    setTxtProgressBar(pb, i_tot)
    i_tot = i_tot + 1
    eta[i_tot,] = r.prop(theta$a.mu[1,], sigma = theta$a.cov[,,1])
    mod = fit.inla(data,eta[i_tot,])
    mlik[i_tot] = mod$mlik
    margs = store.post(mod$dists,margs,i,N_tot)
    delta[i_tot] = N_0*d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1], log = FALSE)
    weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot,]) - d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1]))
  }
  theta$a.mu[2,] = c(sum(eta[1:i_tot,1]*weight[1:i_tot])/sum(weight[1:i_tot]),sum(eta[1:i_tot,2]*weight[1:i_tot])/sum(weight[1:i_tot]))
  theta$a.cov[,,2] = calc.cov(theta$a.mu[2,],weight[1:i_tot],eta[1:i_tot,])
  N_tmp = N_0
  
  # adaptive importance sampling
  for (t in seq(length(N_t))){
    N_tmp = N_tmp + N_t[t]
    for (i in seq(N_t[t])){
      setTxtProgressBar(pb, i_tot)
      i_tot = i_tot + 1
      eta[i_tot,] = r.prop(theta$a.mu[t+1,], sigma = theta$a.cov[,,t+1])
      mod = fit.inla(data,eta[i_tot,])
      mlik[i_tot] = mod$mlik
      margs = store.post(mod$dists,margs,i_tot,N_tot)
      delta[i_tot] = N_0*d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1],log = FALSE) + calc.delta(N_t,eta[i_tot,],theta, t,d.prop)
      weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot,]))/(delta[i_tot]/N_tmp)
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],weight[1:(N_tmp - N_t[t])],N_t = c(N_0,N_t),eta[1:(N_tmp - N_t[t]),],theta,t,mlik[1:(N_tmp - N_t[t])],prior,d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta$a.mu[t+2,] = c(sum(eta[1:i_tot,1]*weight[1:i_tot])/sum(weight[1:i_tot]),sum(eta[1:i_tot,2]*weight[1:i_tot])/sum(weight[1:i_tot]))
    theta$a.cov[,,t+2] = calc.cov(theta$a.mu[t+2,],weight[1:i_tot],eta[1:i_tot,])
  }
  browser()
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
mod = inla.w.is(data = df, init = list(mu = c(0,0), cov = diag(5,2,2)), prior.beta, dq.beta, rq.beta, N_t = rep(20,20))


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

intercept = fit.marginals(mod$weight,mod$margs$intercept)
tau = fit.marginals(mod$weight,mod$margs$tau)
eta.kernel = data.frame(beta1 = mod$eta[,1], beta2 = mod$eta[,2], heat = mod$weight/sum(mod$weight))

ggplot(intercept) + 
  geom_line(aes(x = x, y = y))

ggplot(tau) + 
  geom_line(aes(x = x, y = y))

ggplot(eta.kernel, aes(x = beta1, y = beta2)) + 
  geom_density_2d()
  
  