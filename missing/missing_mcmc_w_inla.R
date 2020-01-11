require(INLA)


dq.x.mis <- function(x, y, sigma = sqrt(0.001), log =TRUE) {
  res <- dnorm(y, mean = x, sd = sigma, log = log)
  
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

rq.x.mis <- function(x, sigma = sqrt(0.001)) {
  rnorm(length(x), mean = x, sd = sigma)
}

mcmc.w.inla <- function(data, init, prior, d.prop, r.prop, fit.inla,
                        n.samples = 100, n.burnin = 5, n.thin = 1){
  eta = matrix(NA, nrow = n.samples, ncol = length(init))
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  eta[1,] = init
  mod.curr = fit.inla(data, eta[1,])
  mlik[1] = mod.curr$mlik
  starttime = Sys.time()
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  i_marg = 0
  N_marg = floor((n.samples - n.burnin)/n.thin)
  times = numeric(N_marg)
  margs = NA
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    eta.new = r.prop(eta[i-1,])
    mod.new = fit.inla(data, eta.new)
    lacc1 = mod.new$mlik + prior(eta.new) + d.prop(eta.new, eta[i-1,])
    lacc2 = mod.curr$mlik + prior(eta[i-1,]) + d.prop(eta[i-1,], eta.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      eta[i,] = eta.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
      
    }else{
      eta[i,] = eta[i-1,]
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i > n.burnin){
      if (((i-1) %% n.thin)==0){
        i_marg = i_marg + 1
        margs = store.post(mod.curr$dists,margs,i_marg,N_marg)
        times[i_marg] =  as.numeric(difftime(Sys.time(),starttime,units = "secs"))
      }
    }
  }
  eta = eta[-seq(n.burnin),]
  eta = eta[seq(from = 1, to = nrow(eta), by=n.thin),]
  return(list(eta = eta,
              margs = lapply(margs, function(x){fit.marginals(rep(1,N_marg),x)}),
              acc.vec = acc.vec,
              mlik = mlik,
              times = times))
}
