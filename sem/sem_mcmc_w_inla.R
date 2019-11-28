

sem.mcmc.w.inla <- function(data, n.samples = 100, n.burnin = 5, n.thin = 1){
  rho = numeric(n.samples)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  mod.curr = fit.inla(data, rho[1])
  mlik[1] = mod.curr$mlik
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    rho.new = rq.rho(rho[i-1])
    mod.new = fit.inla(data,rho.new)
    lacc1 = mod.new$mlik + prior.rho(rho.new) + dq.rho(rho.new, rho[i-1])
    lacc2 = mod.curr$mlik + prior.rho(rho[i-1]) + dq.rho(rho[i-1], rho.new)
    acc = acc = min(1,exp(lacc1 - lacc2))
    if (runif(1)<acc){
      rho[i] = rho.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
    }else{ 
      rho[i] = rho[i-1]
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
  return(list(rho=rho,
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}