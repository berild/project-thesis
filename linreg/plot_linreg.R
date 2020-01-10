library(ggplot2)
library(scales)
library(ggpubr)

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

p1 <-ggplot() +
  #geom_vline(xintercept = ml[1,1]) +  
  #geom_line(data = ml_w_inla_mod$dists$intercept, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.fixed[[1]]), aes(x = x, y = y)) + 
  #geom_line(data = pp_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "IS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x = "",y = "",title=expression(alpha)) + 
  coord_cartesian(xlim = c(min(mcmc_w_inla_mod$margs$intercept$x),max(mcmc_w_inla_mod$margs$intercept$x))) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "linreg_intercept_plot.pdf", plot = p1, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting univariate distribution of the precition; tau

p2 <- ggplot() + 
  #geom_line(data = ml_w_inla_mod$dists$tau, aes(x = x, y = y, color = "ML with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.hyperpar[[1]]), aes(x = x, y = y)) + 
  #geom_line(data = pp_w_inla_mod$margs$tau, aes(x = x, y = y, color = "PP with INLA")) + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(tau)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2
ggsave(filename = "linreg_tau_plot.pdf", plot = p2, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
## Plotting bivariate distribution of model params; betas'

p3 <- ggplot() + 
  geom_point(data = data.frame(x = 2, y = -2), aes(x = x, y = y), shape = 4,size = 3) + 
  geom_text(data = data.frame(x = 2, y = -2), aes(x = x, y = y), label="True Value", vjust=2) + 
  #geom_contour(data = pp_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "PP with INLA")) + 
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) + 
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) + 
  geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) + 
  labs(color = "",x=expression(beta[1]),y=expression(beta[2])) + 
  theme_bw() + 
  theme(legend.position="bottom")
p3
ggsave(filename = "linreg_contour_plot.pdf", plot = p3, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight))
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight/sum(is_w_inla_mod$weight))
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta,mcmc_w_inla_mod$times)

p4 <- ggplot() + 
  geom_hline(yintercept = 10000) + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$ess, aes(x = time, y = ess, color = "MCMC with INLA")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() + 
  labs(color = "",x="Runtime (secs)",y="Effective sample size") + 
  coord_cartesian(xlim=c(min(amis_w_inla_mod$ess$time),max(mcmc_w_inla_mod$ess$time))) + 
  theme_bw() + 
  theme(legend.position="bottom")
p4

ggsave(filename = "linreg_ess_plot.pdf", plot = p4, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


multi_to_uni_kde <- function(multi_kern){
  seq1 = seq(min(multi_kern$x),max(multi_kern$x),length.out = 51)
  seq2 = seq(min(multi_kern$y),max(multi_kern$y),length.out = 51)
  uni_kern1 = data.frame(x = seq1[-51] + 1/2*(seq1[2]-seq1[1]), 
                        y = rep(0,50))
  uni_kern2 = data.frame(x = seq2[-51] + 1/2*(seq2[2]-seq2[1]), 
                        y = rep(0,50))
  int1 = int2 = 0
  for (i in seq(2,51)){
    uni_kern1$y[i-1] = mean(multi_kern$z[(seq1[i-1]<=multi_kern$x)&(seq1[i]>multi_kern$x)])
    int1 = int1 + uni_kern1$y[i-1]*(seq1[i]-seq1[i-1])
    uni_kern2$y[i-1] = mean(multi_kern$z[(seq2[i-1]<=multi_kern$y)&(seq2[i]>multi_kern$y)])
    int2 = int2 + uni_kern2$y[i-1]*(seq2[i]-seq2[i-1])
  }
  uni_kern1$y = uni_kern1$y/int1
  uni_kern2$y = uni_kern2$y/int2
  return(list(kern1 = uni_kern1,
              kern2 = uni_kern2))
}

amis_uni_kerns = multi_to_uni_kde(amis_w_inla_mod$eta_kern)
is_uni_kerns = multi_to_uni_kde(is_w_inla_mod$eta_kern)
mcmc_uni_kerns = multi_to_uni_kde(mcmc_w_inla_mod$eta_kern)

p5 <- ggplot() + 
  geom_line(data = as.data.frame(inla_mod$marginals.fixed$x1), aes(x=x,y=y)) +
  geom_line(data = amis_uni_kerns$kern1, aes(x=x,y=y,color = "AMIS with INLA")) + 
  geom_line(data = is_uni_kerns$kern1, aes(x=x,y=y,color = "IS with INLA")) + 
  geom_line(data = mcmc_uni_kerns$kern1, aes(x=x,y=y, color = "MCMC with INLA")) + 
  coord_cartesian(xlim = c(0.5,3.5)) + 
  labs(color = "", x = "",y="",title=expression(beta[1]))+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggsave(filename = "linreg_beta1_plot.pdf", plot = p5, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p6 <- ggplot() + 
  geom_line(data = as.data.frame(inla_mod$marginals.fixed$x2), aes(x=x,y=y)) + 
  geom_line(data = amis_uni_kerns$kern2, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_uni_kerns$kern2, aes(x=x,y=y, color = "IS with INLA")) + 
  geom_line(data = mcmc_uni_kerns$kern2, aes(x=x,y=y, color = "MCMC with INLA")) + 
  labs(color = "", x = "",y="",title = expression(beta[2]))+
  coord_cartesian(xlim = c(-3.5,-0.5))  +
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
ggsave(filename = "linreg_beta2_plot.pdf", plot = p6, device = NULL, path = "./linreg/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


ptot <- ggarrange(p3, p4, p5, p6, p2, p1, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ptot
ggsave(filename = "linreg_tot.pdf", plot = ptot, device = NULL, path = "./linreg/figures/",
       scale = 1, width = 2*width, height = 3*height, units = "in", dpi=5000)
