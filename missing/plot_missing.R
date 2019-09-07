# INLA within McMC
load(file = "missing-mcmc-w-inla.Rdata")
res = cbind(data.frame(step = seq(nrow(mod$x.mis))),as.data.frame(mod$x.mis))
params = colnames(res[,-c(1)])
res$is_burnin = c(rep(T,5),rep(F,nrow(res)-5))

means = data.frame(key = params, 
                   value = c(sapply(res[res$is_burnin==F,params],mean)))

res = gather(res,key,value,params)
res$key=factor(res$key,levels=params)
res2 = rbind(
  cbind(key = rep("beta0",nrow(mod$beta0)),mod$beta0),
  cbind(key = rep("beta1",nrow(mod$beta1)),mod$beta1),
  cbind(key = rep("beta2",nrow(mod$beta2)),mod$beta2),
  cbind(key = rep("beta3",nrow(mod$beta3)),mod$beta3),
  cbind(key = rep("tau",nrow(mod$tau)),mod$tau)
)

traceplot <- ggplot(res, aes(x = step, y = value, color = is_burnin)) +
  geom_line() +
  facet_wrap(vars(key),scales="free",ncol = 2)

traceplot

distplot <- ggplot(res[res$is_burnin==F,]) +
  geom_density(aes(x = value),fill = "deepskyblue" , alpha = 0.5) +
  geom_vline(data = means,aes(xintercept = value),color = "firebrick") + 
  geom_text(data = means,aes(label = sprintf("%.3f",value), x = value), 
            y = 0.02, size = 6,color = "firebrick") + 
  facet_wrap(vars(key),scales="free",ncol = 2)

distplot

inladist <- ggplot(res2, aes(x = x, y = y)) + 
  geom_line(color = "deepskyblue") + 
  facet_wrap(.~key,scales = "free", ncol = 2)
inladist
