# INLA within McMc
library(ggplot2)
load(file = "./classification/classification-mcmc-w-inla.Rdata")
post.marg = data.frame(x = mod$post.marg[[1]][,1],
                       y = mod$post.marg[[1]][,2], 
                       key = rep(names(mod$post.marg)[1],
                                 length(mod$post.marg[[1]][,1])))
for (i  in seq(2,length(mod$post.marg))){
  post.marg = rbind(post.marg, 
                    data.frame(x = mod$post.marg[[i]][,1],
                               y = mod$post.marg[[i]][,2], 
                               key = rep(names(mod$post.marg)[i],
                                         length(mod$post.marg[[i]][,1]))))
}

ggplot(post.marg, aes(x = x, y = y)) + 
  geom_line()+ 
  facet_wrap(.~key,scales = "free", ncol = 2)

stats = data.frame(key = levels(post.marg$key),
                   mean = sapply(seq(length(mod$post.marg)),
                                 function (X){
                                   mod$post.marg[[X]][which.max(mod$post.marg[[X]][,2]),1]
                                 }))

num.integral <- function(df){
  res = 0
  for(i in seq(2,nrow(df))){
    res = res + (df[i,1]-df[i-1,1])*(df[i,2]+df[i-1,2])*0.5
  }
  res
}

num.integral(mod$post.marg[[4]])

dists =data.frame(x = seq(min(faithful$eruptions),max(faithful$eruptions),0.1),
                  y1 = dnorm(seq(min(faithful$eruptions),max(faithful$eruptions),0.1), 
                             mean = stats$mean[1], 
                             sd = stats$mean[3]),
                  y2 = dnorm(seq(min(faithful$eruptions),max(faithful$eruptions),0.1),
                             mean = stats$mean[2],
                             sd = stats$mean[4]))


ggplot(dists, aes(x = x)) +
  geom_line(aes(y = y1)) 
