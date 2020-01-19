library(ggplot2)
library(scales)
library(ggpubr)

## loading simulation results
load(file = "./linreg/sims/linreg-inla.Rdata")
load(file = "./linreg/sims/linreg-is-w-inla.Rdata")
load(file = "./linreg/sims/linreg-amis-w-inla.Rdata")
load(file = "./linreg/sims/linreg-mcmc-w-inla.Rdata")
width = 7
height = 7

## Plotting univariate distribution of the intercept; alpha

p1 <-ggplot() +
  geom_line(data = data.frame(inla_mod$marginals.fixed[[1]],
                              type = rep("INLA",nrow(inla_mod$marginals.fixed[[1]]))), 
            aes(x = x, y = y,linetype = type)) + 
  geom_vline(xintercept = 3, linetype = "dashed") +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "IS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x = "",y = "",title=expression(alpha),linetype = "") + 
  coord_cartesian(xlim = c(min(mcmc_w_inla_mod$margs$intercept$x),max(mcmc_w_inla_mod$margs$intercept$x))) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "linreg_intercept_plot.pdf", plot = p1, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting univariate distribution of the precition; tau

p2 <- ggplot() + 
  geom_line(data = data.frame(inla_mod$marginals.hyperpar[[1]],
                              type = rep("INLA",nrow(inla_mod$marginals.hyperpar[[1]]))), 
            aes(x = x, y = y,linetype = type)) + 
  geom_vline(xintercept = 1, linetype = "dashed")+
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(tau),linetype = "") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2
ggsave(filename = "linreg_tau_plot.pdf", plot = p2, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting bivariate distribution of model params; betas'

p3 <- ggplot() + 
  geom_point(data = data.frame(x = 2, y = -2), aes(x = x, y = y), shape = 4,size = 3) + 
  geom_text(data = data.frame(x = 2, y = -2), aes(x = x, y = y), label="True Value", vjust=2) + 
  geom_line(data = data.frame(x=rep(1000,10),y = rep(1000,10),type = "INLA"), aes(x=x,y=y,linetype=type))+
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) + 
  labs(color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="") + 
  coord_cartesian(xlim = c(1.2,2.7),ylim=c(-2.7,-1.3))+
  theme_bw() + 
  theme(legend.position="bottom")
p3
ggsave(filename = "linreg_contour_plot.pdf", plot = p3, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta, mcmc_w_inla_mod$times)
save(mcmc_w_inla_mod, file="./linreg/sims/linreg-mcmc-w-inla.Rdata")

p4 <- ggplot() + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$ess, aes(x = time, y = ess, color = "MCMC with INLA")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() + 
  labs(color = "",x="Runtime (sec)",y="Effective sample size") + 
  coord_cartesian(xlim=c(min(amis_w_inla_mod$ess$time),max(mcmc_w_inla_mod$ess$time))) + 
  theme_bw() + 
  theme(legend.position="bottom")
p4

ggsave(filename = "linreg_ess_plot.pdf", plot = p4, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


amis_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
is_kerns = lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_kerns = lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        kernel = "gaussian")[c(1,2)])
})

p5 <- ggplot() + 
  geom_line(data = data.frame(inla_mod$marginals.fixed$x1,
                              type = rep("INLA",nrow(inla_mod$marginals.fixed$x1))), 
            aes(x=x,y=y,linetype = type)) + 
  geom_vline(xintercept = 2, linetype = "dashed") + 
  geom_line(data = amis_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  coord_cartesian(xlim = c(0.5,3.5)) + 
  labs(color = "", x = "",y="",title=expression(beta[1]),linetype = "")+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p5
ggsave(filename = "linreg_beta1_plot.pdf", plot = p5, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p6 <- ggplot() + 
  geom_line(data = data.frame(inla_mod$marginals.fixed$x2, 
                              type = rep("INLA",nrow(inla_mod$marginals.fixed$x2))), 
            aes(x=x,y=y,linetype = type)) + 
  geom_vline(xintercept=-2,linetype = "dashed") + 
  geom_line(data = amis_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "", x = "",y="",title = expression(beta[2]),linetype ="")+
  coord_cartesian(xlim = c(-3.5,-0.5))  +
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p6
ggsave(filename = "linreg_beta2_plot.pdf", plot = p6, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


ptot <- ggarrange(p3, p4, p5, p6, p2, p1, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ptot
ggsave(filename = "linreg_tot.pdf", plot = ptot, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
