library(ggplot2)

## loading simulation results
load(file = "./lasso/sims/lasso-inla.Rdata")
load(file = "./lasso/sims/lasso-ml-w-inla.Rdata")
load(file = "./lasso/sims/lasso-pp-w-inla.Rdata")
load(file = "./lasso/sims/lasso-is-w-inla.Rdata")
load(file = "./lasso/sims/lasso-amis-w-inla.Rdata")
load(file = "./lasso/sims/lasso-mcmc-w-inla.Rdata")


ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA"))+
  labs(color = "") + 
  theme_bw()
