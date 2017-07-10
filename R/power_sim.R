source("../R/util.R")
source("../R/stochastic_model.R")

nsim <- 100

meanlog <- 1.8
sdlog <- 0.25
bU <- 20

nsite <- c(5, 10, 15, 20, 25, 30)

simlist <- array(list(), dim=c(nsim, length(nsite)))

for(i in 1:length(nsite)) {
    n.site <- nsite[i]
    j <- 1
    while(j <= 100) {
        beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
        sim <- stochastic_spatial_discrete_model(beta=beta, bU=bU)
        if (!any(sim$A.count[500:1100,] < 0.1)) {
            simlist[j,i][[1]] <- sim
            j <- j + 1
        }
        cat(j, i, "\n")
        save("simlist", file="powersim.rda")
    }
}

