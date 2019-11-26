## loading simulation results
load(file = "./linreg/sims/linreg-inla.Rdata")
load(file = "./linreg/sims/linreg-ml-w-inla.Rdata")
load(file = "./linreg/sims/linreg-pp-w-inla.Rdata")
load(file = "./linreg/sims/linreg-is-w-inla.Rdata")
load(file = "./linreg/sims/linreg-amis-w-inla.Rdata")
load(file = "./linreg/sims/linreg-mcmc-w-inla.Rdata")


## Plotting univariate distribution of the intercept; alpha

ggplot() +
  geom_vline(xintercept = ml[1,1]) +  
  geom_line(data = ml_w_inla_mod$dists$intercept, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.fixed[[1]]), aes(x = x, y = y, color = "INLA")) + 
  geom_line(data = pp_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "IS with INLA")) + 
  labs(color = "") + 
  theme_bw()

## Plotting univariate distribution of the precition; tau

ggplot() + 
  geom_line(data = ml_w_inla_mod$dists$tau, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.hyperpar[[1]]), aes(x = x, y = y, color = "INLA")) + 
  geom_line(data = pp_w_inla_mod$margs$tau, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  labs(color = "") + 
  theme_bw()

## Plotting bivariate distribution of model params; betas'

ggplot() + 
  geom_point(data = data.frame(x = ml[2,1], y = ml[3,1]), aes(x = x, y = y), shape = 4,size = 3, color = "navyblue") + 
  geom_text(data = data.frame(x = ml[2,1], y = ml[3,1]), aes(x = x, y = y), label="ML", vjust=2) + 
  geom_contour(data = pp_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "PP with INLA")) + 
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA")) + 
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA")) + 
  labs(color = "") + 
  theme_bw()
