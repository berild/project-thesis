library(ggplot2)
## loading simulation results
load(file = "./linreg/sims/linreg-inla.Rdata")
#load(file = "./linreg/sims/linreg-ml-w-inla.Rdata")
#load(file = "./linreg/sims/linreg-pp-w-inla.Rdata")
load(file = "./linreg/sims/linreg-is-w-inla.Rdata")
load(file = "./linreg/sims/linreg-amis-w-inla.Rdata")
load(file = "./linreg/sims/linreg-mcmc-w-inla.Rdata")
width = 5
height = 5

## Plotting univariate distribution of the intercept; alpha

ggplot() +
  #geom_vline(xintercept = ml[1,1]) +  
  #geom_line(data = ml_w_inla_mod$dists$intercept, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.fixed[[1]]), aes(x = x, y = y, color = "INLA")) + 
  #geom_line(data = pp_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "IS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x = "",y = "") + 
  coord_cartesian(xlim = c(min(mcmc_w_inla_mod$margs$intercept$x),max(mcmc_w_inla_mod$margs$intercept$x))) + 
  theme_bw() + 
  theme(legend.position="bottom")
ggsave(filename = "linreg_intercept_plot.pdf", plot = last_plot(), device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting univariate distribution of the precition; tau

ggplot() + 
  #geom_line(data = ml_w_inla_mod$dists$tau, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.hyperpar[[1]]), aes(x = x, y = y, color = "INLA")) + 
  #geom_line(data = pp_w_inla_mod$margs$tau, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x="",y="") + 
  theme_bw() + 
  theme(legend.position="bottom")
ggsave(filename = "linreg_tau_plot.pdf", plot = last_plot(), device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting bivariate distribution of model params; betas'

ggplot() + 
  geom_point(data = data.frame(x = ml[2,1], y = ml[3,1]), aes(x = x, y = y), shape = 4,size = 3) + 
  geom_text(data = data.frame(x = ml[2,1], y = ml[3,1]), aes(x = x, y = y), label="ML", vjust=2) + 
  #geom_contour(data = pp_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "PP with INLA")) + 
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) + 
  labs(color = "",x="",y="") + 
  theme_bw() + 
  theme(legend.position="bottom")
ggsave(filename = "linreg_contour_plot.pdf", plot = last_plot(), device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight))
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight/sum(is_w_inla_mod$weight))
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta)

ggplot() + 
  geom_hline(yintercept = 10000) + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  #geom_line(data = data.frame(time = mcmc_w_inla_mod$times, ess = mcmc_w_inla_mod$ess), aes(x = time, y = ess, linetype = "MCMC with INLA")) + 
  labs(color = "",x="",y="") + 
  theme_bw() + 
  theme(legend.position="bottom")
ggsave(filename = "linreg_ess_plot.pdf", plot = last_plot(), device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

