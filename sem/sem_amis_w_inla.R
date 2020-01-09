library(INLA)
library(spdep)
library(spatialreg)
library(INLABMA)

data(columbus)

lw <- nb2listw(col.gal.nb, style="W")
colsemml <- errorsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
                       quiet=TRUE)
W <- as(as_dgRMatrix_listw(nb2listw(col.gal.nb)), "CsparseMatrix")
columbus$idx<-1:nrow(columbus)
form<- CRIME ~ INC + HOVAL
zero.variance = list(prec=list(initial = 25, fixed=TRUE))
df = columbus

fit.inla <- function(data, rho) {
  res <- sem.inla(form, d = data, W = W, rho = rho,
                  family = "gaussian", impacts = FALSE,
                  control.family = list(hyper = zero.variance),
                  verbose = FALSE)
  
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           INC = res$marginals.fixed[[2]],
                           HOVAL = res$marginals.fixed[[3]],
                           tau = res$marginals.hyperpar[[1]])))
}

dq.rho <- function(x, theta, i_cur,log =TRUE) {
  dnorm(x, mean = theta$a.mu[i_cur], sd = theta$a.cov[i_cur], log = log)
}

rq.rho <- function(theta, i_cur) {
  rnorm(1, mean = theta$a.mu[i_cur], sd = theta$a.cov[i_cur])
}


init.dq.rho <- function(x, log = TRUE) {
  dunif(x, -1.5, 1, log = log)
}

init.rq.rho <- function() {
  runif(1, -1.5, 1)
}

prior.rho <- function(x, log = TRUE) {
  dunif(x, -1.5, 1, log = log)
}

store.post <- function(marg,margs,j,n.prop){
  if (is.na(margs)){
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
    tmp = tmp + N_t[l]*d.prop(x = eta, theta = theta, i_cur = l+1, log = FALSE)
  }
  return(tmp)
}

update.delta.weight <- function(delta,weight,N_t,eta,theta,t,mlik,prior,d.prop){
  i_tmp = 0
  N_tmp = sum(N_t[1:(t+1)])
  for (l in seq(t)){
    for (i in seq(N_t[l])){
      i_tmp = i_tmp + 1
      delta[i_tmp] = delta[i_tmp] + N_t[l]*d.prop(x = eta[i_tmp], theta = theta, i_cur = t+1, log = FALSE)
      weight[i_tmp] = exp(mlik[i_tmp] + prior(eta[i_tmp]))/(delta[i_tmp]/N_tmp)
    }
  }
  return(list(
    delta = delta,
    weight = weight
  ))
}

inla.w.amis <- function(data, prior, d.prop, r.prop, fit.inla, init.r.prop, init.d.prop, N_t = rep(20,20)){
  require(INLA)
  N_0 = round(sum(N_t)/2)
  N_tot = N_0 + sum(N_t)
  mlik = numeric(N_tot)
  eta = rep(NA, N_tot)
  delta = numeric(N_tot)
  weight = numeric(N_tot)
  theta = list(a.mu = rep(NA, length(N_t) + 2),
               a.cov = rep(NA, length(N_t) + 2))
  
  # initialization process 
  i_tot = 0
  pb <- txtProgressBar(min = 0, max = N_tot, style = 3)
  margs = NA
  for (i in seq(N_0)){
    setTxtProgressBar(pb, i_tot)
    i_tot = i_tot + 1
    eta[i_tot] = init.r.prop()
    mod = fit.inla(data,eta[i_tot])
    mlik[i_tot] = mod$mlik
    margs = store.post(mod$dists,margs,i,N_tot)
    delta[i_tot] = N_0*init.d.prop(x = eta[i_tot], log = FALSE)
    weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot]) - init.d.prop(x = eta[i_tot]))
  }
  theta = calc.theta(theta,weight,eta,i_tot,2)
  N_tmp = N_0
  
  # adaptive importance sampling
  for (t in seq(length(N_t))){
    N_tmp = N_tmp + N_t[t]
    for (i in seq(N_t[t])){
      setTxtProgressBar(pb, i_tot)
      i_tot = i_tot + 1
      eta[i_tot] = r.prop(theta,t+1)
      mod = fit.inla(data,eta[i_tot])
      mlik[i_tot] = mod$mlik
      margs = store.post(mod$dists,margs,i_tot,N_tot)
      delta[i_tot] = N_0*init.d.prop(x = eta[i_tot], log = FALSE) + calc.delta(N_t,eta[i_tot],theta,t,d.prop)
      weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot]))/(delta[i_tot]/N_tmp)
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],weight[1:(N_tmp - N_t[t])],N_t = c(N_0,N_t),eta[1:(N_tmp - N_t[t])],theta,t,mlik[1:(N_tmp - N_t[t])],prior,d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta = calc.theta(theta,weight,eta,i_tot,t+2)
  }
  return(list(eta = eta,
              theta = theta,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight))
}

