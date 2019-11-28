require(INLA)
require(mvtnorm)
require(MASS)


lasso.mcmc.w.inla <- function(data, n.beta, stdev.samp, n.samples = 100, n.burnin = 5, n.thin = 1){
  beta = matrix(data = NA,nrow = n.samples, ncol = n.beta)
  colnames(beta) = colnames(data$x)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  beta[1,] = rep(0,n.beta)
  mod.curr = fit.inla(data, beta[1,])
  mlik[1] = mod.curr$mlik
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    beta.new = rq.beta(beta[i-1,],sigma = stdev.samp)
    mod.new = fit.inla(data,beta.new)
    lacc1 = mod.new$mlik + prior.beta(beta.new) + dq.beta(beta.new, beta[i-1,])
    lacc2 = mod.curr$mlik + prior.beta(beta[i-1,]) + dq.beta(beta[i-1,], beta.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      beta[i,] = beta.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
    }else{ 
      beta[i,] = beta[i-1,]
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i == n.burnin){
      post.marg = mod.curr$dists
    }else if(i > n.burnin){
      post.marg = moving.marginals(mod.curr$dists,
                                   post.marg,
                                   i-n.burnin+1)
    }
  }
  return(list(beta = beta, 
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}
