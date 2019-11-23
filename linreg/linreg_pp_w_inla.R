require(INLA)
require(MASS)

pp.w.inla <- function(data, ml, fit.inla ,prior ,len = 10){
  eta = as.matrix(expand.grid(beta1 = seq(from = ml[1,1] - 2*ml[1,2], to = ml[1,1] + 2*ml[1,2], length.out = len),
                              beta2 = seq(from = ml[2,1] - 2*ml[2,2], to = ml[2,1] + 2*ml[2,2], length.out = len)))
  margs = NA
  mlik = numeric(len*len)
  weight = numeric(len*len)
  pb <- txtProgressBar(min = 0, max = nrow(eta), style = 3)
  for (i in seq(nrow(eta))){
    setTxtProgressBar(pb, i)
    mod = fit.inla(data,eta[i,])
    mlik[i] = mod$mlik
    margs = store.post(mod$dists,margs,i,len*len)
    weight[i] = mlik[i] + prior(eta[i,])
  }
  weight = exp(weight - max(weight))
  eta_kern = kde2d.weighted(x = eta[,1], y = eta[,2], w = weight/(sum(weight)), n = 100, lims = c(1,3,-3,-1))
  eta_kern = data.frame(expand.grid(x=eta_kern$x, y=eta_kern$y), z=as.vector(eta_kern$z))
  return(list(eta = eta,
              eta_kern = eta_kern,
              margs = lapply(margs, function(x){fit.marginals(weight,x)}),
              weight = weight))
}
