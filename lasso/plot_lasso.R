# McMC with INLA results
tau = as.data.frame(mod$tau)

ggplot(tau) + 
  geom_line(aes(x= x, y = y))

res = cbind(data.frame(step = seq(nrow(mod$beta)), acc.prob = mod$acc.prob),as.data.frame(mod$beta))
res$is_burnin = c(rep(T,500),rep(F,nrow(res)-500))

ggplot(res) +
  geom_line(aes(x = step, y = acc.prob))

params = colnames(res[,-c(1,2,ncol(res))])

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

# McMC results
res = cbind(data.frame(step = seq(nrow(mod$beta)),acc.prob = mod$acc.prob, tau = mod$tau), as.data.frame(mod$beta))
res$is_burnin = c(rep(T,500),rep(F,nrow(res)-500))
params = colnames(res[,-c(1,2,9)])
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

                   