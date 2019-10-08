library(INLA)
library(ggplot2)
sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}


inla.w.ml <- function(data, ml){
  data$oset = data$x%*%ml
  res = inla(y~1+offset(oset), data = data)
  return(list(intercept = res$marginals.fixed[[1]], 
              tau = res$marginals.hyperpar[[1]]))
}


set.seed(1)
df = sample.linreg()
mod = inla.w.ml(df,ml=c(2,-2))

ggplot(as.data.frame(mod$intercept)) + 
  geom_line(aes(x = x, y = y))


ggplot(as.data.frame(mod$tau)) + 
  geom_line(aes(x = x, y = y))
