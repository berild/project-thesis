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

source("./lasso/lasso_ml_w_inla.R")
is_w_inla_mod = ml.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 6*stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(mod, file = "./lasso/lasso-ml-w-inla.Rdata")


source("./lasso/lasso_pp_w_inla.R")
pp_w_inla_mod = pp.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 6*stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(mod, file = "./lasso/lasso-pp-w-inla.Rdata")


source("./lasso/lasso_is_w_inla.R")
is_w_inla_mod = is.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 6*stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(mod, file = "./lasso/lasso-amis-w-inla.Rdata")


source("./lasso/lasso_amis_w_inla.R")
amis_w_inla_mod = amis.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 6*stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(amis_w_inla_mod, file = "./lasso/lasso-amis-w-inla.Rdata")


source("./lasso/lasso_mcmc_w_inla.Rdata")
mcmc_w_inla_mod = mcmc.w.inla(data = df, init = list(mu = rep(0,n.beta), cov = 6*stdev.samp), prior.beta, dq.beta, rq.beta, fit.inla, N_t = rep(20,20))
save(mod, file = "./lasso/lasso-mcmc-w-inla.Rdata")





















# moving.marginals <- function(marg, post.marg, n){
#   for (i in seq(length(post.marg))){
#     tmp.post.marg = post.marg[[i]]
#     tmp.marg = marg[[i]]
#     step = tmp.post.marg[2,1] - tmp.post.marg[1,1]
#     new.x = tmp.post.marg[,1]
#     new.y = tmp.post.marg[,2]
#     if (max(tmp.marg[,1])>max(new.x)){
#       new.u = seq(from = max(new.x) + step,
#                   to = max(tmp.marg[,1])+ step,
#                   by = step)
#       new.x = c(new.x, new.u)
#       new.y = c(new.y,rep(0,length(new.u)))
#     }
#     if (min(tmp.marg[,1])<min(new.x)){
#       new.l = seq(from = min(new.x) - step,
#                   to = min(tmp.marg[,1]) - step,
#                   by = - step)
#       new.x = c(rev(new.l), new.x)
#       new.y = c(rep(0,length(new.l)),new.y)
#     }
#     tmp.post.marg = data.frame(x = new.x, y = new.y)
#     tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
#     tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
#     post.marg[[i]] = tmp.post.marg
#   }
#   post.marg
# }

