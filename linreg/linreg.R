library(INLA)
library(mvtnorm)
library(MASS)
library(ggplot2)


sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

set.seed(1)
df = sample.linreg()