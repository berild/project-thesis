library(MASS)#For 'geyser' dataset
library(MCMCpack)#For dirichlet distribution
library(INLA)
library(parallel)
options(mc.cores = 2)


yy <- faithful$eruptions

n <- length(yy)