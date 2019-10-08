
sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 3 + 2*x1 -2*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}


moving.marginals <- function(marg, post.marg, n){
  for (i in seq(length(post.marg))){
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
    tmp.post.marg[,2] = (tmp + (n-1)*tmp.post.marg[,2])/n
    post.marg[[i]] = tmp.post.marg
  }
  post.marg
}

fit.inla = function(data, ccd.params){
  data$oset = data$x%*%ccd.params
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]])))
}

create.ccd <-function(ml,sd,n = 10){
  params = length(ml)
  ccd.seq = matrix(NA, nrow = params, ncol = n)
  for (i in seq(params)){
    ccd.seq[i,] = seq(ml[i]-0.75*sd[i],ml[i]+0.75*sd[i], length.out = n)
  }
  return(ccd.seq)
}

inla.w.ccd <- function(data,ml,sd,n.samples = 10){
  ccd_seq = create.ccd(ml,sd,n.samples)
  mlik = numeric(n.samples)
  for (i in seq(n.samples)){
    mod = fit.inla(data,ccd_seq[,i])
    mlik[i] = mod$mlik
    if (i == 1){
      post.marg = mod$dists
    }else{
      post.marg = moving.marginals(mod$dists,post.marg,i)
    }
  }
  return(list(beta = ccd_seq,
              post.marg = post.marg,
              mlik = mlik))
  
}

set.seed(1)
df = sample.linreg()
mod = inla.w.ccd(df,ml=c(2,-2),as.numeric(summary(lm(y~x,data=df))$coefficients[-1,2]),n.samples = 100)



ggplot(as.data.frame(mod$post.marg$intercept)) + 
  geom_line(aes(x = x, y = y))


ggplot(as.data.frame(mod$post.marg$tau)) + 
  geom_line(aes(x = x, y = y))
