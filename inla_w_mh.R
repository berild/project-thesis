moving.marginals <- function(marg, post.marg, n){
  for (i in seq(length(post.marg))){
    tmp.post.marg = post.marg[[i]]
    tmp.marg = marg[[i]]
    step = tmp.post.marg[2,1] - tmp.post.marg[1,1]
    new.x = tmp.post.marg[,1]
    new.y = tmp.post.marg[,2]
    if (max(tmp.marg[,1])>max(new.x)){
      new.u = seq(from = max(new.x) + step,
                  to = max(tmp.marg[,1])+ step,
                  by = step)
      new.x = c(new.x, new.u)
      new.y = c(new.y,rep(0,length(new.u)))
    }
    if (min(tmp.marg[,1])<min(new.x)){
      new.l = seq(from = min(new.x) - step,
                  to = min(tmp.marg[,1]) - step,
                  by = - step)
      new.x = c(rev(new.l), new.x)
      new.y = c(rep(0,length(new.l)),new.y)
    }
    tmp.post.marg = data.frame(x = new.x, y = new.y)
    tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
    tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
    post.marg[[i]] = tmp.post.marg
  }
  post.marg
}

inla.w.mcmc <- function(data, init_fixed, prior.fixed, rq.fixed, dq.fixed, 
                         fit.inla, n.samples = 100, n.burnin = 5, n.thin = 1){
  require(INLA)
  fixed = matrix(data = NA,nrow = n.samples, ncol = ncol(data)-1)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  colnames(fixed) = colnames(fixed[,-1])
  fixed[1,] = init_fixed
  mod.curr = fit.inla(data,fixed[1,])
  mlik[1] = mod.curr$mlik
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2, n.samples)){
    setTxtProgressBar(pb, i)
    fixed.new = rq.fixed(fixed[i-1,])
    mod.new = fit.inla(data, fixed.new)
    lacc1 = mod.new$mlik + prior.fixed(fixed.new) + dq.fixed(fixed.new, fixed[i-1,])
    lacc2 = mod.curr$mlik + prior.fixed(beta[i-1,]) + dq.fixed(fixed[i-1,], fixed.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      fixed[i,] = fixed.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
      
    }else{
      fixed[i,] = fixed[i-1,]
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
  return(list(fixed = fixed,
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}


