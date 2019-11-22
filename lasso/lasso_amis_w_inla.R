library(INLA)
library(ggplot2)
library(ISLR)
library(glmnet)
library(smoothmest)
library(mvtnorm)


data(Hitters)

Hitters <- na.omit(Hitters)

x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5]
y <- Hitters$Salary


y <- scale(y)
x <- scale(x)
df <- list(y = y, x = x)
n.beta <- ncol(df$x)

# finding inverse of the precision
stdev.samp <- .25 * solve(t(x)%*%x)

fit.inla <- function(data, b) {
  data$oset = data$x %*% b
  res = inla(y ~ -1 + offset(oset), data = data)
  res = inla.rerun(res)
  return(list(mlik = res$mlik[[1]],  
              dists = list(tau = res$marginals.hyperpar[[1]])))
}

prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))
  
  if(!log) { res <- exp(res) }
  
  return(res)
}

dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
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

kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n) 
  gy <- seq(lims[3], lims[4], length = n) 
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] 
  ay <- outer(gy, y, "-")/h[2]
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) 
  return(list(x = gx, y = gy, z = z))
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

inla.w.amis <- function(data, init, prior, d.prop, r.prop, fit.inla,N_t = rep(20,20)){
  require(INLA)
  N_0 = round(sum(N_t)/2)
  N_tot = N_0 + sum(N_t)
  mlik = numeric(N_tot)
  eta = matrix(NA, ncol = length(init$mu), nrow = N_tot)
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
      delta[i_tot] = N_0*d.prop(y = eta[i_tot,], x = theta$a.mu[1,], sigma = theta$a.cov[,,1],log = FALSE) + calc.delta(N_t,eta[i_tot,],theta, t,d.prop)
      weight[i_tot] = exp(mlik[i_tot] + prior(eta[i_tot,]))/(delta[i_tot]/N_tmp)
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],weight[1:(N_tmp - N_t[t])],N_t = c(N_0,N_t),eta[1:(N_tmp - N_t[t]),],theta,t,mlik[1:(N_tmp - N_t[t])],prior,d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta = calc.theta(theta,weight,eta,i_tot,t+2)
  }
  return(list(eta = eta,
              theta = theta,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight))
}

set.seed(123)
mod = inla.w.amis(data = df, init = list(mu = rep(0,n.beta), cov = stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(mod, file = "./lasso/lasso_amis_w_inla.Rdata")

mod_inla = inla(y~x, data = df)
res_inla = list(intercept = data.frame(mod_inla$marginals.fixed[[1]]),
                cov1 = data.frame(mod_inla$marginals.fixed[[2]]),
                cov2 = data.frame(mod_inla$marginals.fixed[[3]]),
                tau = data.frame(mod_inla$marginals.hyperpar[[1]]))

ggplot() + 
  geom_line(data = mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA"))

ggplot() + 
  geom_contour(data=mod$eta_kern,aes(x=x,y=y,z=z,color = "AMIS with INLA"))


