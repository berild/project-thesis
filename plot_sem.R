rho = data.frame(rho = mod$rho)
rho2 =data.frame(rho = rho[-seq(50000),])
intercept = as.data.frame(mod$intercept)
INC = as.data.frame(mod$INC)
HOVAL = as.data.frame(mod$HOVAL)
tau = as.data.frame(mod$tau)

ggplot(rho2, aes(x = rho)) + 
  geom_density()


ggplot(INC, aes(x = x, y = y)) + 
  geom_line()

ggplot(HOVAL, aes(x = x, y = y)) + 
  geom_line()

ggplot(tau, aes(x = x, y = y)) + 
  geom_line()

ggplot(intercept, aes(x = x, y = y)) + 
  geom_line()
