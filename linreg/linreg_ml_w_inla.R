require(INLA)
fit.inla.ml <- function(data,ml){
  data$oset = data$x%*%ml[-1,1]
  res = inla(y~offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = data.frame(res$marginals.fixed[[1]]), 
                           tau = data.frame(res$marginals.hyperpar[[1]]))))
}
