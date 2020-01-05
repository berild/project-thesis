require(INLA)
require(MASS)
require(mvtnorm)
require(parallel)

rq.beta <- function(x=c(0,0), sigma = diag(5,2,2)) {
  rmvnorm(1,mean=x,sigma = sigma)
}

dq.beta <- function(y, x, sigma = diag(5,2,2), log =TRUE) {
  dmvnorm(y,mean = x, sigma = sigma,log = log)
}

par.is <- function(x, data, theta, t, prior, d.prop, r.prop, fit.inla){
  eta = r.prop(theta$a.mu[t,], sigma = theta$a.cov[,,t])
  mod = fit.inla(data = data ,beta = eta)
  weight = mod$mlik + prior(eta) - d.prop(y = eta, x = theta$a.mu[t,], sigma = theta$a.cov[,,t])
  return(list(mlik = mod$mlik, dists = mod$dists, stats = mod$stats, eta = eta, weight = weight, times = Sys.time()))
}


is.w.inla <- function(data, init, prior, d.prop, r.prop, N_0 = 200, N = 400){
  eta = matrix(NA, ncol = length(init$mu), nrow = N_0)
  weight = numeric(N_0)
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), 2))) 
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  starttime = Sys.time()
  times = numeric(N)
  pb <- txtProgressBar(min = 0, max = N+N_0, style = 3)
  is.list = mclapply(seq(N_0), function(x){
    par.is(x, data, theta, 1, prior, d.prop,r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = detectCores())
  for (i in seq(length(is.list))){
    setTxtProgressBar(pb, i)
    eta[i,] = is.list[[i]]$eta
    weight[i] = is.list[[i]]$weight
  }
  weight = exp(weight - max(weight))
  theta = calc.theta(theta,weight,eta,nrow(eta),2)
  margs = NA
  stats = NA
  eta = matrix(NA, ncol = length(init$mu), nrow = N)
  weight = numeric(N)
  mlik = numeric(N)
  is.list = mclapply(seq(N), function(x){
    par.is(x, data, theta, 2, prior, d.prop,r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = detectCores())
  for (i in seq(length(is.list))){
    setTxtProgressBar(pb, i+N_0)
    eta[i,] = is.list[[i]]$eta
    margs = store.post(is.list[[i]]$dists,margs,i,N)
    stats = store.stats(is.list[[i]]$stats,stats,i,N)
    weight[i] = is.list[[i]]$weight
    mlik[i] = is.list[[i]]$mlik
    times[i] = as.numeric(is.list[[i]]$times-starttime)
  }
  weight = exp(weight - max(weight))
  eta_kern = kde2d.weighted(x = eta[,1], y = eta[,2], w = weight/(sum(weight)), n = 100, lims = c(1,3,-3,-1))
  eta_kern = data.frame(expand.grid(x=eta_kern$x, y=eta_kern$y), z=as.vector(eta_kern$z))
  return(list(eta = eta,
              theta = theta,
              eta_kern = eta_kern,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight,
              times = times))
}

