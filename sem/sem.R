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

dq.rho <- function(x, y, sigma = .15, log =TRUE) {
  dnorm(y, mean = x, sd = sigma, log = log)
}

rq.rho <- function(x, sigma = .15) {
  rnorm(1, mean = x, sd = sigma)
}

prior.rho <- function(x, log = TRUE) {
  dunif(x, -1.5, 1, log = log)
}



dist.init <- function(x,n.step = 100){
  lx = min(x)
  ux = max(x)
  step = (ux - lx)/(n.step-1)
  seq(from = lx - step, to = ux + step, by = (ux - lx)/(n.step-1))
}


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

mod <- sem.mcmc.w.inla(columbus,n.samples = 100000, n.burnin = 500)
save(mod, file = "sem-mcmc-w-inla.Rdata")

