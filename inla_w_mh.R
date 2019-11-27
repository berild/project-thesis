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

moving.bma <- function(marg, post.marg, n){
  require(INLA)
  for (i in seq(length(marg))){
    if(length(post.marg)!=length(marg)){
      return(marg)
    }else{
      tmp.post.marg = post.marg[[i]]
      tmp.marg = marg[[i]]
      min.x = min(c(tmp.post.marg[1,1],tmp.marg[1,1]))
      max.x = max(c(tmp.post.marg[length(tmp.post.marg[,1]),1], tmp.marg[length(tmp.marg[,1]),1]))
      new.x = seq(min.x,max.x,length.out = 100)
      new.y = inla.dmarginal(new.x, tmp.post.marg, log = FALSE)
      # new.x = tmp.post.marg[,1]
      # new.y = tmp.post.marg[,2]
      # if ((max(tmp.marg[,1])>max(new.x)) & (min(tmp.marg[,1])<(min(new.x)))){
      #   
      # }else if (max(tmp.marg[,1])>max(new.x)){
      #   new.u = seq(from = max(new.x) + step,
      #               to = max(tmp.marg[,1])+ step,
      #               by = step)
      #   new.x = c(new.x, new.u)
      #   new.y = c(new.y,rep(0,length(new.u)))
      # }else if (min(tmp.marg[,1])<min(new.x)){
      #   new.l = seq(from = min(new.x) - step,
      #               to = min(tmp.marg[,1]) - step,
      #               by = - step)
      #   new.x = c(rev(new.l), new.x)
      #   new.y = c(rep(0,length(new.l)),new.y)
      # }
      tmp.post.marg = data.frame(x = new.x, y = new.y)
      tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
      tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
      post.marg[[i]] = tmp.post.marg
    }
  }
  return(post.marg)
}

inla.w.mcmc <- function(data, init_fixed, prior.fixed, rq.fixed, dq.fixed, 
                         fit.inla, n.samples = 100, n.burnin = 5, n.thin = 1){
  require(INLA)
  fixed = matrix(data = NA,nrow = n.samples, ncol = ncol(data$x))
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  colnames(fixed) = colnames(data$x[,-1])
  fixed[1,] = init_fixed
  mod.curr = fit.inla(data,fixed[1,])
  mlik[1] = mod.curr$mlik
  post.marg = NA
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2, n.samples)){
    setTxtProgressBar(pb, i)
    fixed.new = rq.fixed(fixed[i-1,])
    mod.new = fit.inla(data, fixed.new)
    lacc1 = mod.new$mlik + prior.fixed(fixed.new) + dq.fixed(fixed.new, fixed[i-1,])
    lacc2 = mod.curr$mlik + prior.fixed(fixed[i-1,]) + dq.fixed(fixed[i-1,], fixed.new)
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
    if (i > n.burnin)
      post.marg = moving.bma(mod.curr$dists,post.marg,i-n.burnin)
  }
  return(list(fixed = fixed,
              post.marg = post.marg,
              acc.vec = acc.vec,
              mlik = mlik))
}


