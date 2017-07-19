source("../R/util.R")
source("../R/stochastic_model.R")

nsim <- 100

meanlog <- seq(0.2, 3, by=0.4)
sdlog <- 0.23
bU <- 20.6

nsite <- 20

simlist <- array(list(), dim=c(nsim, length(meanlog)))

for(i in 1:length(meanlog)) {
    j <- 1
    while(j <= 100) {
        beta <- rlnorm(nsite, meanlog=meanlog[i], sdlog=sdlog)
        sim <- stochastic_spatial_discrete_model(beta=beta, bU=bU,n.site=nsite)
        if (!any(sim$A.count[500:1100,] < 0.1)) {
            simlist[j,i][[1]] <- sim
            j <- j + 1
        }
        cat(j, i, "\n")
        save("simlist", file="betasim.rda")
    }
}
