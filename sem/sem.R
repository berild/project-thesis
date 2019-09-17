library(INLA)
library(spdep)
library(spatialreg)
library(INLABMA)

data(columbus)

lw <- nb2listw(col.gal.nb, style="W")
colsemml <- errorsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
                       quiet=FALSE)
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
              intercept = res$marginals.fixed[[1]],
              INC = res$marginals.fixed[[2]],
              HOVAL = res$marginals.fixed[[3]],
              tau = res$marginals.hyperpar[[1]]))
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


moving.marginals <- function(new.marginal, avg.marginal, n){
  step = avg.marginal$x[2] - avg.marginal$x[1]
  new.x = avg.marginal$x
  new.y = avg.marginal$y
  if (max(new.marginal[,"x"])>max(new.x)){
    new.u = seq(from = max(new.x) + step,
                to = max(new.marginal[,"x"])+ step,
                by = step)
    new.x = c(new.x, new.u)
    new.y = c(new.y,rep(0,length(new.u)))
  }
  if (min(new.marginal[,"x"])<min(new.x)){
    new.l = seq(from = min(new.x) - step,
                to = min(new.marginal[,"x"]) - step,
                by = - step)
    new.x = c(rev(new.l), new.x)
    new.y = c(rep(0,length(new.l)),new.y)
  }
  avg.marginal = data.frame(x = new.x, y = new.y)
  tmp = inla.dmarginal(avg.marginal$x, new.marginal, log = FALSE)
  avg.marginal$y = (tmp + (n-1)*avg.marginal$y)/n
  avg.marginal
}

sem.mcmc.w.inla <- function(data, n.samples = 100, n.burnin = 5, n.thin = 1){
  rho = numeric(n.samples)
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  mod.curr = fit.inla(data, rho[1])
  mlik[1] = mod.curr$mlik
  intercept =data.frame(x = rep(0,102), y = rep(0,102))
  INC = data.frame(x = rep(0,102), y = rep(0,102))
  HOVAL = data.frame(x = rep(0,102), y = rep(0,102))
  tau = data.frame(x = rep(0,102), y = rep(0,102))
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    rho.new = rq.rho(rho[i-1])
    mod.new = fit.inla(data,rho[i])
    lacc1 = mod.new$mlik + prior.rho(rho.new) + dq.rho(rho.new, rho[i-1])
    lacc2 = mod.curr$mlik + prior.rho(rho[i-1]) + dq.rho(rho[i-1], rho.new)
    acc = acc = min(1,exp(lacc1 - lacc2))
    paste0("acceptance prob", acc)
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
      intercept$x = dist.init(x = mod.curr$intercept[,"x"])
      INC$x = dist.init(x = mod.curr$INC[,"x"])
      HOVAL$x = dist.init(x = mod.curr$HOVAL[,"x"])
      tau$x = dist.init(x = mod.curr$tau[,"x"])
    }else if(i > n.burnin){
      intercept = moving.marginals(new.marginal = mod.curr$intercept,
                               avg.marginal = intercept,
                               n = i-n.burnin)
      INC = moving.marginals(new.marginal = mod.curr$INC,
                               avg.marginal = INC,
                               n = i-n.burnin)
      HOVAL = moving.marginals(new.marginal = mod.curr$HOVAL,
                               avg.marginal = HOVAL,
                               n = i-n.burnin)
      tau = moving.marginals(new.marginal = mod.curr$tau,
                             avg.marginal = tau,
                             n = i-n.burnin)
    }
  }
  return(list(rho=rho,
              intercept =intercept, 
              INC = INC, 
              HOVAL = HOVAL,
              tau = tau,
              acc.vec = acc.vec))
}

mod <- sem.mcmc.w.inla(columbus,n.samples = 100000, n.burnin = 500)
save(mod, file = "sem-mcmc-w-inla.Rdata")

