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

sem.mcmc.w.inla <- function(data){
  N = 100000
  burnin = 500
  rho = c(0)
  mod1 = fit.inla(data, rho[1])
  intercept = mod1$intercept*0
  INC = mod1$INC*0
  HOVAL = mod1$HOVAL*0
  tau = mod1$tau*0
  acc.prob = c()
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  for (i in seq(2,N)){
    setTxtProgressBar(pb, i)
    rho = c(rho,rq.rho(rho[i-1]))
    mod2 = fit.inla(data,rho[i])
    acc.prob = c(acc.prob,
                 mod2$mlik + 
                   prior.rho(rho[i]) +
                   dq.rho(rho[i],rho[i-1]) - 
                   mod1$mlik -
                   prior.rho(rho[i-1]) - 
                   dq.rho(rho[i-1], rho[i]))
    if (log(runif(1))>acc.prob[i-1]){
      rho[i] = rho[i-1]
      if (i > burnin){
        intercept = intercept + mod1$intercept
        INC = INC + mod1$INC
        HOVAL = HOVAL + mod1$HOVAL
        tau = tau + mod1$tau
      }
    }else if (i > burnin){
      intercept = intercept + mod2$intercept
      INC = INC + mod2$INC
      HOVAL = HOVAL + mod2$HOVAL
      tau = tau + mod2$tau
    }
  }
  return(list(rho=rho,
              intercept = intercept/(N-burnin), 
              INC = INC/(N-burnin), 
              HOVAL = HOVAL/(N-burnin),
              tau = tau/(N-burnin),
              acc.prob = sapply(exp(acc.prob),min,1)))
}

mod <- sem.mcmc.w.inla(columbus)
save(mod, file = "sem-mcmc-w-inla.Rdata")

