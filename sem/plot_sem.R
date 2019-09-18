# McMC with INLA results
load(file = "./sem/sem-mcmc-w-inla.Rdata")
rho = data.frame(x = mod$rho)

post.marg <- data.frame(x = c(mod$intercept$x,
                              mod$INC$x,
                              mod$HOVAL$x,
                              mod$tau$x),
                        y = c(mod$intercept$y,
                              mod$INC$y,
                              mod$HOVAL$y,
                              mod$tau$y),
                        key = c(rep("intercept",length(mod$intercept$x)),
                                rep("INC",length(mod$INC$x)),
                                rep("HOVAL",length(mod$HOVAL$x)),
                                rep("tau",length(mod$tau$x))))


ggplot(rho, aes(x = x)) + 
  geom_density(color = "deepskyblue", fill = "deepskyblue", alpha = 0.3)

ggplot(post.marg, aes(x = x, y = y)) + 
  geom_line(color = "deepskyblue") + 
  facet_wrap(.~key,scales = "free", ncol = 2)

mean(mod$acc.vec)
