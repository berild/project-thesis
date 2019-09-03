
ggplot(as.data.frame(mod$tau),aes(x = x, y = y)) +
  geom_line()
as.data.frame(mod$beta)
cbind(data.frame(step = seq(1,nrow(mod$beta))),as.data.frame(mod$beta))
ggplot(cbind(data.frame(step = seq(1,nrow(mod$beta))),as.data.frame(mod$beta))) +
  geom_line(aes(x = step, y = V2, color = "beta_1")) +
  labs(color = "")

beta = as.data.frame(mod$beta)
colnames(beta) = c("beta_1","beta_2","beta_3","beta_4","beta_5")
params = colnames(beta)
beta2 = beta[-seq(5000),]
ggplot(beta2,aes(x = beta_4)) +
  geom_density(fill = "deepskyblue", alpha = 0.5)

beta %>%
  gather(key,value, params) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "deepskyblue" , alpha = 0.5) +
  facet_wrap(vars(key),scales="free",ncol = 1)