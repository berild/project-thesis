library(ggplot2)
library(scales)
library(ggpubr)

## loading simulation results
load(file = "./sem/sims/sem-is-w-inla.Rdata")
load(file = "./sem/sims/sem-amis-w-inla.Rdata")
load(file = "./sem/sims/sem-mcmc-w-inla.Rdata")


marginal.stat <- function(method){
  require(INLA)
  stats = lapply(method, function(y){
    m1 = inla.emarginal(function(x) x,y)
    m2 = inla.emarginal(function(x) x^2,y)
    c(m1,sqrt(m2-m1^2))
  })
  return(stats)
}

amis_w_inla_mod$stats = marginal.stat(amis_w_inla_mod$margs)
is_w_inla_mod$stats = marginal.stat(is_w_inla_mod$margs)
mcmc_w_inla_mod$stats = marginal.stat(mcmc_w_inla_mod$margs)

width = 7
height = 7

p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,color="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title="Intercept") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(25,95))
p1
ggsave(filename = "sem_intercept.pdf", plot = p1, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p2 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$INC, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$INC, aes(x=x,y=y,color="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$INC, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title="INC") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(-2.8,1))
p2
ggsave(filename = "sem_INC.pdf", plot = p2, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$HOVAL, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$HOVAL, aes(x=x,y=y,color="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$HOVAL, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title="HOVAL") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(-1,0.4))
p3
ggsave(filename = "sem_HOVAL.pdf", plot = p3, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,color="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(tau)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(0.004,0.02))
p4
ggsave(filename = "sem_tau.pdf", plot = p4, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

amis_kern = as.data.frame(density(x = amis_w_inla_mod$eta,
                                  weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight), 
                                  kernel = "gaussian", from = -1, to = 2)[c(1,2)])

is_kern = as.data.frame(density(x = is_w_inla_mod$eta,
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight), 
                        kernel = "gaussian")[c(1,2)])

mcmc_kerns = lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        kernel = "gaussian")[c(1,2)])
})

p5 <-ggplot() + 
  geom_line(data=amis_kern, aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kern, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(rho)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim = c(0,1))
p5
ggsave(filename = "sem_rho_plot.pdf", plot = p5, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

source("./sem/sem_general_functions.R")

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta,mcmc_w_inla_mod$times)

p6 <- ggplot() + 
  #geom_hline(yintercept = 10000) + 
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  #geom_line(data = mcmc_w_inla_mod$ess, aes(x = time, y = ess, color = "MCMC with INLA")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() + 
  labs(color = "",x="Runtime (sec)",y="Effective sample size") + 
  #coord_cartesian(xlim=c(min(amis_w_inla_mod$ess$time),max(mcmc_w_inla_mod$ess$time))) + 
  theme_bw() + 
  theme(legend.position="bottom")
p6
ggsave(filename = "sem_ess_plot.pdf", plot = p6, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

ptot <- ggarrange(p1, p2, p3, p4, p5, p6, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ptot
ggsave(filename = "sem_tot.pdf", plot = ptot, device = NULL, path = "./sem/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

