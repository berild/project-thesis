library(ggplot2)
library(scales)
library(ggpubr)

library(ISLR)
library(glmnet)
library(coda)

data(Hitters)
Hitters <- na.omit(Hitters)
#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)


#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = 0.073,intercept = F)
lasso.coef


## loading simulation results
load(file = "./lasso/sims/lasso-is-w-inla.Rdata")
load(file = "./lasso/sims/lasso-amis-w-inla.Rdata")
load(file = "./lasso/sims/lasso-mcmc-w-inla.Rdata")

width = 7
height = 7

p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(tau)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "lasso_tau_plot.pdf", plot = p1, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


amis_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})
is_kerns = lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        bw = 0.08,
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})
mcmc_kerns = lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})


p2 <-ggplot() + 
  geom_vline(xintercept = lasso.coef[2],linetype = "dashed") + 
  geom_line(data = amis_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="AtBat") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2
ggsave(filename = "lasso_b1_plot.pdf", plot = p2, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p3 <-ggplot() + 
  geom_vline(xintercept = lasso.coef[3],linetype = "dashed") + 
  geom_line(data = amis_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="Hits") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p3
ggsave(filename = "lasso_b2_plot.pdf", plot = p3, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p4 <-ggplot() + 
  geom_vline(xintercept = lasso.coef[4],linetype = "dashed") + 
  geom_line(data = amis_kerns[[3]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[3]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[3]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="HmRun") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p4
ggsave(filename = "lasso_b3_plot.pdf", plot = p4, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p5 <-ggplot() + 
  geom_vline(xintercept = lasso.coef[5],linetype = "dashed") + 
  geom_line(data = amis_kerns[[4]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[4]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[4]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="Runs") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p5
ggsave(filename = "lasso_b4_plot.pdf", plot = p5, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p6 <-ggplot() + 
  geom_vline(xintercept = lasso.coef[6],linetype = "dashed") + 
  geom_line(data = amis_kerns[[5]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[5]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[5]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="RBI") + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p6
ggsave(filename = "lasso_b5_plot.pdf", plot = p6, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight))
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight/sum(is_w_inla_mod$weight))

p7 <- ggplot() + 
  #geom_hline(yintercept = 10000) + 
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
p7
ggsave(filename = "lasso_ess_plot.pdf", plot = p7, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

ptot <- ggarrange(p2, p3, p4, p5, p6, p7, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ptot
ggsave(filename = "lasso_tot.pdf", plot = ptot, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
