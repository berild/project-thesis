library(MASS)#For 'geyser' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
INLA:::inla.dynload.workaround()
library(parallel)
options(mc.cores = 2)


yy <- faithful$eruptions

n <- length(yy)
n.grp <- 2
grp <- rep(2, n)
grp[order(yy)[1:floor(n/3)]] <- 1


fit.inla.internal <- function(yy, grp) {
  y = matrix(NA, ncol = n.grp, nrow = n)
  for(i in 1:n.grp) {
    idx = which(grp == i)
    y[idx, i] = yy[idx]
  }
  x = y
  x[!is.na(x)] = 1
  d = list(y = y, x = x)
  m1 = inla(y ~ -1 + x, data = d,
             family = rep("gaussian", n.grp),
             control.fixed = list(mean = list(x1 = 2, x2 = 4.5), prec = 1)
  )
  return(list(means = m1$summary.fixed[, "mean"],
              precs = m1$summary.hyperpar[, "mean"], 
              mlik = m1$mlik[[1]],
              dists = list(mu1 = as.data.frame(m1$marginals.fixed[[1]]),
                           mu2 = as.data.frame(m1$marginals.fixed[[2]]),
                           tau1 = as.data.frame(m1$marginals.hyperpar[[1]]),
                           tau2 = as.data.frame(m1$marginals.hyperpar[[2]]))))
}

get.probs <- function(z) {
  probs = rep(0, n.grp)
  tab = table(z)
  probs[as.integer(names(tab))] = tab/sum(tab)
  return(probs)
}

dq.z <- function(x, y, log = TRUE) {
  means = x$m$means
  precs = x$m$means
  ww = get.probs(x$z)
  z.probs = sapply(1:n, function (X) {
    aux = ww * dnorm(yy[X], means, sqrt(1/precs))
    (aux/sum(aux))[y$z[X]]
  })
  if(log) {
    return(sum(log(z.probs)))
  } else {
    return(prod(z.probs))
  }
}

rq.z <- function(z) {
  means = z$m$means
  precs = z$m$precs
  probs = get.probs(z$z)
  z.sim = sapply(1:n, function (X) {
    aux = probs * dnorm(yy[X], means, sqrt(1/precs))
    sample(1:n.grp, 1, prob = aux/sum(aux))
  })
  z.model = fit.inla.internal(yy, z.sim)
  z.new = list(z = z.sim, m = z.model)
  return(z.new)
}


prior.z <- function(z, log = TRUE) {
  res <- log(0.5) * length(z$z)
  if(log) {
    return(res)
  }
  else {
    return(exp(res))
  }
}


moving.marginals <- function(marg, mlik, post.marg, n){
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
    tmp.post.marg[,2] = (tmp*exp(mlik) + (n-1)*tmp.post.marg[,2])/n
    post.marg[[i]] = tmp.post.marg
  }
  post.marg
}


classification.mcmc.w.inla <- function(data, grp, n.samples = 100, n.burnin = 5, n.thin = 1){
  grp.curr = list(z = grp, m = fit.inla.internal(yy, grp))
  grp.curr = rq.z(grp.curr)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  z = matrix(NA,nrow=n.samples,ncol=length(grp))
  z[1,] = grp.curr$z
  mlik[1] = grp.curr$m$mlik
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2, n.samples)){
    setTxtProgressBar(pb, i)
    grp.new = rq.z(grp.curr)
    lacc1 = grp.new$m$mlik + prior.z(grp.new) + dq.z(grp.new, grp.curr)
    lacc2 = grp.curr$m$mlik + prior.z(grp.curr) + dq.z(grp.curr, grp.new)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      z[i,] = grp.new$z
      grp.curr = grp.new
      mlik[i] = grp.new$m$mlik
      acc.vec[i] = T
      
    }else{
      z[i,] = z[i-1,]
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if (n.burnin == i){
      post.marg = grp.curr$m$dists*grp.curr$m$mlik
    }else if (i > n.burnin){
      post.marg = moving.marginals(grp.curr$m$dists,
                                   grp.curr$m$mlik,
                                   post.marg,n = i-n.burnin + 1)
    }
  }
  return(list(z = z,
              z.probs = apply(z,2,get.probs),
              acc.vec = acc.vec,
              post.marg = post.marg,
              mlik = mlik))
}

mod = classification.mcmc.w.inla(faithful$eruptions,grp,n.samples = 100000, n.burnin = 500, n.thin = 1)
save(mod, file = "./classification/classification-mcmc-w-inla.Rdata")