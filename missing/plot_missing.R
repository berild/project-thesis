# INLA within McMC
res = cbind(data.frame(step = seq(nrow(mod$x.mis)), acc.prob = mod$acc.prob),as.data.frame(mod$x.mis))
params = colnames(res[,-c(1,2)])
res$is_burnin = c(rep(T,500),rep(F,nrow(res)-500))

means = data.frame(key = params, 
                   value = c(sapply(res[,params],mean)))
res = gather(res,key,value,params)
res$key=factor(res$key,levels=params)
res2 = rbind(
  cbind(key = rep("beta_0",nrow(mod$beta_0)),as.data.frame(beta_0)),
  cbind(key = rep("beta_1",nrow(mod$beta_1)),as.data.frame(beta_1)),
  cbind(key = rep("beta_2",nrow(mod$beta_2)),as.data.frame(beta_2)),
  cbind(key = rep("beta_3",nrow(mod$beta_3)),as.data.frame(beta_3)),
  cbind(key = rep("tau",nrow(mod$tau)),as.data.frame(tau))
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