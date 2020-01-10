library(mice) # data
library(INLA)
library(mvtnorm)
library(parallel)
library(coda)
library(MASS)

data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)

df = list(d.mis = d.mis, idx.mis = idx.mis)

init = list(mu = rep(mean(df$d.mis$bmi, na.rm = TRUE),n.mis),
            cov = diag(2*var(df$d.mis$bmi, na.rm = TRUE),n.mis,n.mis))

names(init$mu) = sprintf("Observation_%d",df$idx.mis)
colnames(init$cov) = sprintf("Observation_%d",df$idx.mis)

source("./missing/missing_general_function.R")

#source("./missing/missing_amis_w_inla.R")
#amis_w_inla_mod <- amis.w.inla(data = df, init = init, prior.x.mis, 
#                               dq.x.mis, rq.x.mis, fit.inla, 
#                               N_t = seq(25,50,1)*10, N_0 = 250)
#save(amis_w_inla_mod, file = "./missing/sims/missing-amis-w-inla.Rdata")
#
#source("./missing/missing_is_w_inla.R")
#is_w_inla_mod <- is.w.inla(data = df, init = init, prior.x.mis, 
#                           dq.x.mis, rq.x.mis,fit.inla, N_0 = 800, N = 10000)
#save(is_w_inla_mod, file = "./missing/sims/missing-is-w-inla.Rdata")

source("./missing/missing_mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = init$mu,
                               prior.x.mis, dq.x.mis, rq.x.mis, fit.inla,
                               n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./missing/sims/missing-mcmc-w-inla.Rdata")


