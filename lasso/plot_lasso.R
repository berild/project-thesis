library(ggplot2)

## loading simulation results
#load(file = "./lasso/sims/lasso-inla.Rdata")
#load(file = "./lasso/sims/lasso-ml-w-inla.Rdata")
#load(file = "./lasso/sims/lasso-pp-w-inla.Rdata")
load(file = "./lasso/sims/lasso-is-w-inla.Rdata")
load(file = "./lasso/sims/lasso-amis-w-inla.Rdata")
load(file = "./lasso/sims/lasso-mcmc-w-inla.Rdata")

width = 5
height = 5

p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) + 
  #geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(tau)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "lasso_tau_plot.pdf", plot = p2, device = NULL, path = "./lasso/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


multi_to_uni_kde <- function(eta,weight=NA,interval=NA,binwidth=50){
  if (anyNA(weight)){
    weight = 1
  }else{
    weight = weight/sum(weight)
  }
  if (anyNA(interval)){
    xseq = seq(min(eta),max(eta),length.out = binwidth+1) 
  }else{
    xseq = seq(interval[1],interval[2],length.out = binwidth+1)
  }
  uni_kern = data.frame(x = xseq[-(binwidth+1)] + 1/2*(xseq[2]-xseq[1]),
                        y = rep(0,binwidth))
  int_kern = 0
  for (j in seq(2,binwidth+1)){
    uni_kern$y[j-1] = sum(weight[(xseq[j-1]<=eta)&(xseq[j]>eta)])
    int_kern = int_kern + uni_kern$y[j-1]*(xseq[j]-xseq[j-1])
  }
  uni_kern$y = uni_kern$y/int_kern
  uni_spline = smooth.spline(uni_kern)
  uni_kern = as.data.frame(predict(uni_spline,x=seq(xseq[1],xseq[binwidth+1],length.out = 4*binwidth)))
  return(uni_kern)
}

amis_kern1 = multi_to_uni_kde(amis_w_inla_mod$eta[,1],amis_w_inla_mod$weight,c(-0.5,0.5),binwidth=15)
amis_kern2 = multi_to_uni_kde(amis_w_inla_mod$eta[,2],amis_w_inla_mod$weight,c(-0.3,0.6),binwidth=15)
amis_kern3 = multi_to_uni_kde(amis_w_inla_mod$eta[,3],amis_w_inla_mod$weight,c(-0.3,0.4),binwidth=15)
amis_kern4 = multi_to_uni_kde(amis_w_inla_mod$eta[,4],amis_w_inla_mod$weight,c(-0.2,0.4),binwidth=15)
amis_kern5 = multi_to_uni_kde(amis_w_inla_mod$eta[,5],amis_w_inla_mod$weight,c(-0.2,0.6),binwidth=15)

p2 <-ggplot() + 
  geom_line(data=amis_kern1, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data=is_kern1, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kern1, aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(beta[1])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2


p3 <-ggplot() + 
  geom_line(data=amis_kern2, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data=is_kern1, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kern1, aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(beta[2])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p3

p4 <-ggplot() + 
  geom_line(data=amis_kern3, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data=is_kern1, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kern1, aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(beta[3])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p4

p5 <-ggplot() + 
  geom_line(data=amis_kern4, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data=is_kern1, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kern1, aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(beta[4])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p5

p6 <-ggplot() + 
  geom_line(data=amis_kern5, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data=is_kern1, aes(x=x,y=y,color="IS with INLA")) + 
  #geom_line(data=mcmc_kern1, aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=expression(beta[5])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p6

