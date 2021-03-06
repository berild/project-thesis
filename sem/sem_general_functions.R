require(INLA)
require(INLABMA)
require(coda)
require(spdep)
require(spatialreg)
require(mvtnorm)
require(MASS)

prior.rho <- function(x, log = TRUE) {
  dunif(x, -1.5, 1, log = log)
}

fit.inla <- function(data, eta) {
  res <- sem.inla(form, d = data, W = W, rho = eta,
                  family = "gaussian", impacts = FALSE,
                  control.family = list(hyper = zero.variance),
                  verbose = FALSE)
  
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           INC = res$marginals.fixed[[2]],
                           HOVAL = res$marginals.fixed[[3]],
                           tau = res$marginals.hyperpar[[1]])))
}

calc.theta <- function(theta,weight,eta,i_tot,i_cur){
    theta$a.mu[i_cur] = sum(eta[1:i_tot]*weight[1:i_tot])/sum(weight[1:i_tot])
    theta$a.cov[i_cur] = sum(weight[1:i_tot]*(eta[1:i_tot]-theta$a.mu[i_cur])*
                               (eta[1:i_tot]-theta$a.mu[i_cur]))/(sum(weight[1:i_tot]))
  return(theta)
}

store.post <- function(marg,margs,j,n.prop){
  if (anyNA(margs)){
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

running.ESS <- function(eta, times, ws = NA, norm = TRUE,step = 100){
  if (anyNA(ws)){
    require(coda)
    ess = sapply(seq(2,length(eta)),function(x){
      effectiveSize(eta[1:x])
    })
    times = times[-1]
  }else{
    if (norm){
      ws = ws/sum(ws)
    }
    ess = sapply(seq(length(ws)),function(x){
      sum(ws[1:x])^2/(sum(ws[1:x]^2))
    })
    rm.ess = !is.na(ess)
    times = times[rm.ess]
    ess = ess[rm.ess]
  }
  ess.df = data.frame(time = c(times[1],times[rev(seq(length(times),100,-step))]),
                      ess = c(ess[1],ess[rev(seq(length(ess),100,-step))]))
  return(ess.df)
}

fit.marginals <- function(ws,margs,len = 400){
  ws = ws/sum(ws)
  xmin <- quantile(apply(margs[[1]],1,function(X){min(X)}),0.25)
  xmax <- quantile(apply(margs[[1]],1,function(X){max(X)}),0.75)
  xx <- seq(xmin, xmax, len = len)
  marg = numeric(len)
  for (i in seq(nrow(margs[[1]]))){
    marg = marg + ws[i]*inla.dmarginal(xx, list(x = margs[[1]][i,], y = margs[[2]][i,]))
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
