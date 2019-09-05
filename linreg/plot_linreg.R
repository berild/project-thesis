# McMC with INLA results
load(file = "./linreg/linreg.Rdata")
tau = as.data.frame(mod$tau)
alfa = as.data.frame(mod$alfa)
res = cbind(data.frame(step = seq(nrow(mod$beta))),as.data.frame(mod$beta))
res$is_burnin = c(rep(T,100),rep(F,nrow(res)-100))

params = colnames(res[,-c(1,ncol(res))])

means = data.frame(key = params, 
                   value = c(sapply(res[,params],mean)))

res = gather(res,key,value,params)

traceplot <- ggplot(res, aes(x = step, y = value, color = is_burnin)) +
  geom_line() +
  facet_wrap(vars(key),scales="free",ncol = 2)

traceplot

distplot <- ggplot(res[res$is_burnin==F,]) +
  geom_density(aes(x = value),fill = "deepskyblue" , alpha = 0.5) +
  geom_vline(data = means,aes(xintercept = value),color = "firebrick") + 
  geom_text(data = means,aes(label = sprintf("%.3f",value), x = value), 
            y = 0.2, size = 6,color = "firebrick") + 
  facet_wrap(vars(key),scales="free",ncol = 2)

distplot

# INLA restults
load(file = "./linreg/linreg_INLA.Rdata")
beta_0 = as.data.frame(mod_inla$marginals.fixed[[1]]) 
beta_1 = as.data.frame(mod_inla$marginals.fixed[[2]])
beta_2 = as.data.frame(mod_inla$marginals.fixed[[3]])
tau = as.data.frame(mod_inla$marginals.hyperpar[[1]])

res = rbind(
  cbind(key = rep("beta_0",nrow(beta_0)),beta_0),
  cbind(key = rep("beta_1",nrow(beta_1)),beta_1),
  cbind(key = rep("beta_2",nrow(beta_2)),beta_2),
  cbind(key = rep("tau",nrow(tau)),tau)
)

rm(mod_inla)

inla_dens <- ggplot(res, aes(x = x, y = y)) + 
  geom_line(color = "deepskyblue") + 
  facet_wrap(.~key, scales = "free", ncol = 2) + 
  theme_classic()
inla_dens
