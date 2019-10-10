rq.beta <- function(x, sigma = .75) {
  rnorm(length(x), mean = x, sd = sigma)
}

dq.beta <- function(x, sigma = .75, log =TRUE) {
  sum(dnorm(x, mean = 0, sd = sigma, log = log))
}

prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}


fit.inla = function(data, beta){
  data$oset = data$x%*%beta
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]])))
}


moving.marginals <- function(marg, post.marg = NA, weight){
  for (i in seq(length(marg))){
    if(length(post.marg)!=length(marg)){
      marg[[i]][,2] = marg[[i]][,2]*weight
      if (i == length(marg)){
        post.marg = marg
      }
    }else{
      tmp.post.marg = post.marg[[i]]
      tmp.marg = marg[[i]]
      step = tmp.post.marg[2,1] - tmp.post.marg[1,1]
      new.x = tmp.post.marg[,1]
      new.y = tmp.post.marg[,2]
      if (max(tmp.marg[,1])>max(new.x)){
        new.u = seq(from = max(new.x) + step,
                    to = max(tmp.marg[,1])+ step,
                    by = step)
        new.x = c(new.x, new.u)
        new.y = c(new.y,rep(0,length(new.u)))
      }
      if (min(tmp.marg[,1])<min(new.x)){
        new.l = seq(from = min(new.x) - step,
                    to = min(tmp.marg[,1]) - step,
                    by = - step)
        new.x = c(rev(new.l), new.x)
        new.y = c(rep(0,length(new.l)),new.y)
      }
      tmp.post.marg = data.frame(x = new.x, y = new.y)
      tmp = inla.dmarginal(tmp.post.marg[,1], tmp.marg, log = FALSE)
      tmp.post.marg[,2] = tmp*weight + tmp.post.marg[,2]
      post.marg[[i]] = tmp.post.marg
    }
  }
  post.marg
}


# implement adaptive mean and standard diviation
# E(\beta) = \sum \beta_i*w_i
# something for standard diviation
# draw samples from this distribution
# have a very vague proposal before



inla.w.is <- function(data, prior, d.prop, n.prop){
  prop.space = matrix(rep(seq(-10,10,1),n.prop), 
                      nrow = n.prop)
  mlik = numeric(ncol(prop.space))
  weight = numeric(ncol(prop.space))
  post.marg = NA
  for (i in seq(ncol(prop.space))){
    mod = fit.inla(data,prop.space[,i])
    mlik[i] = mod$mlik
    weight[i] = exp(mlik[i] + prior(prop.space[,i]) - d.prop(prop.space[,i]))
    post.marg = moving.marginals(mod$dists,post.marg,weight[i])
  }
  
  return(post.marg)
}

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
mod = inla.w.is(df,prior.beta, dq.beta, n.prop = 2)


ggplot(as.data.frame(mod$intercept)) + 
  geom_line(aes(x = x, y = y))


ggplot(as.data.frame(mod$tau)) + 
  geom_line(aes(x = x, y = y))