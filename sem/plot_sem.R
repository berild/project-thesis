# McMC with INLA results
load(file = "./sem/sem-mcmc-w-inla.Rdata")
rho = data.frame(x = mod$rho)


post.marg <- data.frame(x = c(mod$post.marg$intercept$x[mod$post.marg$intercept$y>10^(-3)],
                              mod$post.marg$INC$x[mod$post.marg$INC$y>10^(-3)],
                              mod$post.marg$HOVAL$x[mod$post.marg$HOVAL$y>10^(-3)],
                              mod$post.marg$tau$x[mod$post.marg$tau$y>10^(-3)]),
                        y = c(mod$post.marg$intercept$y[mod$post.marg$intercept$y>10^(-3)],
                              mod$post.marg$INC$y[mod$post.marg$INC$y>10^(-3)],
                              mod$post.marg$HOVAL$y[mod$post.marg$HOVAL$y>10^(-3)],
                              mod$post.marg$tau$y[mod$post.marg$tau$y>10^(-3)]),
                        key = c(rep("intercept",length(mod$post.marg$intercept$x[mod$post.marg$intercept$y>10^(-3)])),
                                rep("INC",length(mod$post.marg$INC$x[mod$post.marg$INC$y>10^(-3)])),
                                rep("HOVAL",length(mod$post.marg$HOVAL$x[mod$post.marg$HOVAL$y>10^(-3)])),
                                rep("tau",length(mod$post.marg$tau$x[mod$post.marg$tau$y>10^(-3)]))))


ggplot(rho, aes(x = x)) + 
  geom_density(color = "deepskyblue", fill = "deepskyblue", alpha = 0.3)

ggplot(post.marg, aes(x = x, y = y)) + 
  geom_line(color = "deepskyblue") + 
  facet_wrap(.~key,scales = "free", ncol = 2)

mean(mod$acc.vec)

stats = data.frame(key = levels(post.marg$key),
                   mean = sapply(seq(length(mod$post.marg)),
                                 function (X){
                                   mod$post.marg[[X]][which.max(mod$post.marg[[X]][,2]),1]
                                 }))
