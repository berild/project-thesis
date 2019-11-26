require(INLA)
require(MASS)
require(mvtnorm)

rq.beta <- function(x=c(0,0), sigma = diag(5,2,2)) {
  rmvnorm(1,mean=x,sigma = sigma)
}

dq.beta <- function(y, x, sigma = diag(5,2,2), log =TRUE) {
  dmvnorm(y,mean = x, sigma = sigma,log = log)
}

calc.theta <- function(theta,weight,eta,i_tot,i_cur){
  for (i in seq(ncol(eta))){
    theta$a.mu[i_cur,i] = sum(eta[1:i_tot,i]*weight[1:i_tot])/sum(weight[1:i_tot])
  }
  for (i in seq(ncol(eta))){
    for (j in seq(i,ncol(eta))){
      theta$a.cov[i,j,i_cur] = theta$a.cov[j,i,i_cur] = sum(weight[1:i_tot]*(eta[1:i_tot,i]-theta$a.mu[i_cur,i])*(eta[1:i_tot,j]-theta$a.mu[i_cur,j]))/(sum(weight[1:i_tot]))
    }
  }
  return(theta)
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
