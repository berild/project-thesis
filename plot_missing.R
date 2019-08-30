
x.mis = as.data.frame(mod$x.mis)
colnames(x.mis) = c("Observation_1","Observation_3","Observation_4",
                    "Observation_6","Observation_10","Observation_11",
                    "Observation_12","Observation_16","Observation_21")
params = colnames(x.mis)
x.mis2 = x.mis[-seq(500),]

intercept = as.data.frame(mod$alfa)
beta = as.data.frame(mod$beta)
tau = as.data.frame(mod$tau)

ggplot(intercept,aes(x = x,y = y)) + 
  geom_line() + 
  labs(title = "Intercept")

ggplot(beta, aes(x = x, y = y)) + 
  geom_line() + 
  labs(title = "Beta")

ggplot(tau, aes(x = x, y = y)) + 
  geom_line() + 
  labs(title = "tau")

x.mis2 %>%
  gather(key,value, params) %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "deepskyblue" , alpha = 0.5) +
  facet_wrap(vars(key),scales="free",ncol = 1)

ggplot(x.mis2) + 
  geom_density(aes(x = Observation_1),fill = "deepskyblue", alpha = 0.5)
