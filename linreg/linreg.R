# loading required libraries
library(INLA)
library(mvtnorm)
library(MASS)
library(parallel)
library(coda)

# function for generating samples
sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

# sampling dataset
set.seed(1)
df = sample.linreg()

# finding inla solution
inla_mod = inla(y~x, data=df) 
save(inla_mod, file = "./linreg/sims/linreg-inla.Rdata")

# finding maximum likelihood
ml = summary(lm(y~x, data = df))$coefficients[,1:2]

source("./linreg/linreg_general_functions.R")

# fitting inla conditioned on samples from is
source("./linreg/linreg_is_w_inla.R")
is_w_inla_mod <- is.w.inla(data = df, init = list(mu = c(0,0), cov = diag(5,2,2)), 
                           prior.beta, dq.beta, rq.beta, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./linreg/sims/linreg-is-w-inla.Rdata")

# fitting inla conditioned on samples from amis
source("./linreg/linreg_amis_w_inla.R")
amis_w_inla_mod <- amis.w.inla(data = df, init = list(mu = c(0,0), cov = diag(5,2,2)), 
                               prior.beta, dq.beta, rq.beta, fit.inla, N_t = seq(25,50,1)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./linreg/sims/linreg-amis-w-inla.Rdata")

# fitting inla conditioned on samples from mh
source("./linreg/linreg_mcmc_w_inla.R")
mcmc_w_inla_mod <- mcmc.w.inla(data = df, init = c(0,0), 
                               prior.beta, dq.beta, rq.beta, fit.inla, 
                               n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./linreg/sims/linreg-mcmc-w-inla.Rdata")
