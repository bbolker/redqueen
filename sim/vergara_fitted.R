source("../R/util.R")
source("../R/stochastic_model.R")

nsim <- 500

meanlog <- 1.82
sdlog <- 0.16
bU <- 20.8
epsilon.site <- 0.007
V <- 0.85
bI <- bU * (1-V)
n.site <- 50

simlist <- vector("list", nsim)

j <- 1

while(j <= 100) {
    beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
    sim <- stochastic_spatial_discrete_model(beta=beta, bU=bU, bI=bI,
                                             epsilon.site=epsilon.site,
                                             n.site=n.site)
    if (!any(sim$A.count[500:1100,] < 0.1)) {
        simlist[[j]] <- sim
        j <- j + 1
    }
    cat(j, i, "\n")
    save("simlist", file="vergara_fitted.rda")
}


