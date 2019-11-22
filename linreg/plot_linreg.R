# AMIS with INLA
load(file = "./linreg/linreg_amis_w_inla.Rdata")
mod_amis = mod
load(file = "./linreg/linreg_mcmc_w_inla.Rdata")
mod_mcmc = mod
mod_inla = inla(y~x, data = df)
mod_inla = list(intercept = data.frame(mod_inla$marginals.fixed[[1]]),
                cov1 = data.frame(mod_inla$marginals.fixed[[2]]),
                cov2 = data.frame(mod_inla$marginals.fixed[[3]]),
                tau = data.frame(mod_inla$marginals.hyperpar[[1]]))

ggplot() + 
  geom_line(data = mod$margs$intercept,aes(x = x, y = y , color = "AMIS with INLA")) + 
  geom_line(data = res_inla$intercept, aes(x=x,y=y, color = "INLA"))

ggplot() + 
  geom_line(data = mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = res_inla$tau, aes(x = x, y = y, color = "INLA"))

ggplot() + 
  geom_contour(data=mod$eta_kern,aes(x=x,y=y,z=z,color = "AMIS with INLA"))
