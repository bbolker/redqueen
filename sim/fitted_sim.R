source("../R/util.R")
source("../R/stochastic_model.R")

L <- load("../data/SMC_summary.rda")

simfun <- function(nsim=30, 
                   beta.meanlog,
                   beta.sdlog,
                   epsilon.site,
                   n.genotype,
                   V, 
                   c_b,
                   n.site=30) {
    simlist <- vector("list", nsim)
    j <- 1
    while(j <= nsim) {
        beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
        bU <- 20
        bI <- bU * (1-V)
        sim <- stochastic_spatial_discrete_lim_model(
            beta=beta,
            n.site=n.site,
            n.genotype=round(n.genotype),
            epsilon.site=epsilon.site,
            bU=bU,
            bI=bI,
            c_b=c_b
        )
        ## avoid fixation
        if(!(all(sim$S.count[1000:1100] < 10) || all(sim$A.count[1000:1100] < 10))) {
            simlist[[j]] <- sim
            j <- j + 1
        }
        cat(j, "\n")
    }
    simlist
}

vergara_sim <- do.call(simfun, SMC_summary$parlist[4,3:8])
save("vergara_sim", file="fitted_sim.rda")

dagan_sim <- do.call(simfun, SMC_summary$parlist[3,3:8])
save("vergara_sim", "dagan_sim", file="fitted_sim.rda")

