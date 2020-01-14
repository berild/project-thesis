require(INLA)
require(spdep)
require(spatialreg)
require(INLABMA)
require(parallel)
require(mvtnorm)

dq.rho <- function(y, x, sigma = init$cov,log =TRUE) {
  dmvt(y,sigma = matrix(sigma), df=3, delta = x, type = "shifted")
  #dnorm(y, mean = x, sd = sigma, log = log)
}

rq.rho <- function(x, sigma = init$cov) {
  as.vector(rmvt(1,sigma = matrix(sigma), df=3, delta = x, type = "shifted"))
  #rnorm(1, mean = x, sd = sigma)
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
  eta = matrix(eta)
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

par.amis <- function(x, data, theta, t, N_0, N_t, N_tmp,
                     prior, d.prop, r.prop, fit.inla){
  eta = r.prop(theta$a.mu[t+1,], sigma = theta$a.cov[,,t+1])
  mod = fit.inla(data = data ,eta = eta)
  if (t==0){
    delta = N_0*d.prop(y = eta, x = theta$a.mu[1,], sigma = theta$a.cov[,,1], log = FALSE)
    weight = exp(mod$mlik + prior(eta) - d.prop(y = eta, x = theta$a.mu[1,], sigma = theta$a.cov[,,1]))
  }else{
    delta = N_0*d.prop(y = eta, x = theta$a.mu[1,], sigma = theta$a.cov[,,1],log = FALSE) + calc.delta(N_t,eta,theta, t, d.prop)
    weight = exp(mod$mlik + prior(eta))/(delta/N_tmp)
  }
  return(list(mlik = mod$mlik, dists = mod$dists, eta = eta, delta = delta, weight = weight, times = Sys.time()))
}


amis.w.inla <- function(data, init, prior, d.prop, r.prop, fit.inla, N_t = rep(20,20), N_0 = NA){
  if (anyNA(N_0)){
    N_0 = round(sum(N_t)/2) 
  }
  if (detectCores()>10){
    ncores = 10
  }else{
    ncores = detectCores()
  }
  N_tot = N_0 + sum(N_t)
  mlik = numeric(N_tot)
  eta = matrix(NA, ncol = length(init$mu), nrow = N_tot)
  delta = numeric(N_tot)
  weight = numeric(N_tot)
  times = numeric(N_tot)
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = length(N_t) + 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), length(N_t) +2)))
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  # initialization process 
  i_tot = 0
  pb <- txtProgressBar(min = 0, max = N_tot, style = 3)
  margs = NA
  starttime = Sys.time()
  N_tmp = N_0
  t = 0
  amis.list = mclapply(seq(N_0),function(x){
    par.amis(x, data, theta, t, N_0, 
             N_t, N_tmp, prior, d.prop,
             r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = ncores)
  for (ele in amis.list){
    setTxtProgressBar(pb, i_tot)
    i_tot = i_tot + 1
    margs = store.post(ele$dists,margs,i_tot,N_tot)
    eta[i_tot,] = ele$eta
    mlik[i_tot] = ele$mlik
    delta[i_tot] = ele$delta
    weight[i_tot] = ele$weight
    times[i_tot] = as.numeric(difftime(ele$times,starttime,units = "secs"))
    
  }
  theta = calc.theta(theta,weight,eta,i_tot,2)
  # adaptive importance sampling
  for (t in seq(length(N_t))){
    N_tmp = N_tmp + N_t[t]
    amis.list = mclapply(seq(N_t[t]),function(x){
      par.amis(x, data, theta, t, N_0, 
               N_t, N_tmp, prior, d.prop,
               r.prop, fit.inla)
    },mc.set.seed = TRUE, mc.cores = ncores)
    for (ele in amis.list){
      setTxtProgressBar(pb, i_tot)
      i_tot = i_tot + 1
      margs = store.post(ele$dists,margs,i_tot,N_tot)
      eta[i_tot,] = ele$eta
      mlik[i_tot] = ele$mlik
      delta[i_tot] = ele$delta
      weight[i_tot] = ele$weight
      times[i_tot] = as.numeric(difftime(ele$times,starttime,units = "secs"))
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],weight[1:(N_tmp - N_t[t])],N_t = c(N_0,N_t),eta[1:(N_tmp - N_t[t]),],theta,t,mlik[1:(N_tmp - N_t[t])],prior,d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta = calc.theta(theta,weight,eta,i_tot,t+2)
  }
  return(list(eta = eta,
              theta = theta,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight,
              times = times))
}
