library(INLA)
library(spdep)
library(spData)
library(spatialreg)
library(INLABMA)
library(parallel)
library(mvtnorm)
library(MASS)
library(coda)

data(columbus)

lw <- nb2listw(col.gal.nb, style="W")
colsemml <- errorsarlm(CRIME ~ INC + HOVAL, data=columbus, lw, method="eigen", 
                       quiet=TRUE)
W <- as(as_dgRMatrix_listw(nb2listw(col.gal.nb)), "CsparseMatrix")
columbus$idx<-1:nrow(columbus)
form<- CRIME ~ INC + HOVAL
zero.variance = list(prec=list(initial = 25, fixed=TRUE))
df = columbus
init = list(mu = 0, cov = 2)

source("./sem/sem_general_functions.R")

#source("./sem/sem_amis_w_inla.R")
#amis_w_inla_mod <- amis.w.inla(data = df, init = init, prior.rho, 
#                               dq.rho, rq.rho, fit.inla, 
#                               N_t = seq(25,50,1)*10, N_0 = 250)
#save(amis_w_inla_mod, file = "./sem/sims/sem-amis-w-inla.Rdata")

#source("./sem/sem_is_w_inla.R")
#is_w_inla_mod <- is.w.inla(data = df, init = init, prior.rho, 
#                           dq.rho, rq.rho,fit.inla, N_0 = 800, N = 10000)
#save(is_w_inla_mod, file = "./sem/sims/sem-is-w-inla.Rdata")

source("./sem/sem_mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init$mu,
                               prior.rho, dq.rho, rq.rho, fit.inla,
                               n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./sem/sims/sem-mcmc-w-inla.Rdata")
