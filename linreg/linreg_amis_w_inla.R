require(mvtnorm)
require(INLA)
require(MASS)

rq.beta <- function(x=c(0,0), sigma = diag(5,2,2)) {
  rmvnorm(1,mean=x,sigma = sigma)
}

dq.beta <- function(y, x, sigma = diag(5,2,2), log =TRUE) {
  dmvnorm(y,mean = x, sigma = sigma,log = log)
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

amis.w.inla <- function(data, init, prior, d.prop, r.prop, fit.inla,N_t = rep(20,20)){
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
  stats = NA
  for (i in seq(N_0)){
    setTxtProgressBar(pb, i_tot)
    i_tot = i_tot + 1
    eta[i_tot,] = r.prop(theta$a.mu[1,], sigma = theta$a.cov[,,1])
    mod = fit.inla(data,eta[i_tot,])
    mlik[i_tot] = mod$mlik
    margs = store.post(mod$dists,margs,i_tot,N_tot)
    stats = store.stats(mod$stats,stats,i_tot,N_tot)
    delta[i_tot] = N_0*d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1], log = FALSE)
    weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot,]) - d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1]))
  }
  theta = calc.theta(theta,weight,eta,i_tot,2)
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
      stats = store.stats(mod$stats,stats,i_tot,N_tot)
      delta[i_tot] = N_0*d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1],log = FALSE) + calc.delta(N_t,eta[i_tot,],theta, t,d.prop)
      weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot,]))/(delta[i_tot]/N_tmp)
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],weight[1:(N_tmp - N_t[t])],N_t = c(N_0,N_t),eta[1:(N_tmp - N_t[t]),],theta,t,mlik[1:(N_tmp - N_t[t])],prior,d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta = calc.theta(theta,weight,eta,i_tot,t+2)
  }
  stats = calc.stats(stats,weight)
  eta_kern = kde2d.weighted(x = eta[,1], y = eta[,2], w = weight/(sum(weight)), n = 100, lims = c(1,3,-3,-1))
  eta_kern = data.frame(expand.grid(x=eta_kern$x, y=eta_kern$y), z=as.vector(eta_kern$z))
  return(list(eta = eta,
              theta = theta,
              eta_kern = eta_kern,
              stats = stats,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight))
}
