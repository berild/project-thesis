# loading required packages
library(INLA)
library(ISLR)
library(glmnet)
library(smoothmest)
library(mvtnorm)
library(parallel)
library(coda)

data(Hitters)
summary(Hitters)

#Check NA's and fix
sum(is.na(Hitters$Salary))
Hitters <- na.omit(Hitters)

#
# The Lasso
#

#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)
n.beta <- ncol(df$x)

# ml estimates
ml = summary(lm(y~x, data = df))$coefficients[2:6,1:2]

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)


#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid)
plot(lasso.mod)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1)
plot(cv.out)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])
mean((lasso.pred - y[test])^2)

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)

#Check estimated coefficients
lasso.coef
lasso.coef[lasso.coef != 0]

#Fitted values
lasso.fitted <- predict(out, s = bestlam, newx = x)
# importing dataset

# finding inverse of the precision
stdev.samp <- .25 * solve(t(x)%*%x)

source("./lasso/lasso_general_function.R")

# fitting inla conditioned on samples from is
source("./lasso/lasso_is_w_inla.R")
is_w_inla_mod = is.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 4*stdev.samp), 
                          prior.beta, dq.beta, rq.beta, fit.inla, N_0 = 800, N = 10000)
save(is_w_inla_mod, file = "./lasso/sims/lasso-is-w-inla.Rdata")

# fitting inla conditioned on samples from amis
source("./lasso/lasso_amis_w_inla.R")
amis_w_inla_mod = amis.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 4*stdev.samp), 
                              prior.beta, dq.beta, rq.beta, fit.inla, N_t = seq(25,50,1)*10, N_0 = 250)
save(amis_w_inla_mod, file = "./lasso/sims/lasso-amis-w-inla.Rdata")

# fitting inla conditioned on samples from mcmc
source("./lasso/lasso_mcmc_w_inla.R")
mcmc_w_inla_mod = mcmc.w.inla(data = df, init = rep(0,ncol(df$x)), 
                              prior.beta, dq.beta, rq.beta, fit.inla, 
                              n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_w_inla_mod, file = "./lasso/sims/lasso-mcmc-w-inla.Rdata")

