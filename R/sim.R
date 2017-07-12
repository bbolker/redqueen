source("util.R")
source("model.R")

n.site <- 20
#q <- (2*(1:n.site)-1)/(2*n.site)

beta_vec <- rlnorm(n.site, meanlog = 0.5, sdlog=0.3)

res <- spatial_discrete_model(beta=beta_vec, n.site=n.site, seed=113, tmax=1100, tburnin=500)

S <- res$S.count[1000,]
A <- res$A.count[1000,]
I <- (res$SI.count+res$AI.count)[1000,]

cor.test((I/(S+A)), (S/(S+A)))

## playing around...
source("stochastic_model.R")

res_spatial <- spatial_discrete_model(beta=5, n.site=1, tmax=1100, tburnin=500, epsilon.site = 0,seed=101)
res <- discrete_model(beta=5, tmax=1100, tburnin = 500, seed=101)
res2 <- stochastic_spatial_discrete_model(beta=5, n.site = 1, seed=103)
res3 <- stochastic_discrete_model(beta=5, seed=103)

plot(res_spatial$S.count, type="l")
lines(res$S.count, col=2)
lines(res2$S.count, col=3)
lines(res3$S.count, col=4)
