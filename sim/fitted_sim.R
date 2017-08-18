source("../R/three_loci_util.R")
source("../R/three_loci_stochastic_model.R")

L <- load("../data/SMC_summary.rda")

simfun <- function(nsim=50, 
                   beta.meanlog,
                   beta.sdlog,
                   epsilon.site,
                   n.genotype,
                   V, 
                   n.site=50) {
    simlist <- vector("list", nsim)
    j <- 1
    while(j <= nsim) {
        beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
        bU <- 20
        bI <- bU * (1-V)
        sim <- three_loci_stochastic_spatial_discrete_lim_model(
            beta=beta,
            n.site=n.site,
            n.genotype=round(n.genotype),
            epsilon.site=epsilon.site,
            bU=bU,
            bI=bI
        )
        
        if (!any(sim$A.count[1000:1100,] < 0.1)) {
            simlist[[j]] <- sim
            j <- j + 1
        }
        cat(j, "\n")
    }
    simlist
}

vergara_sim <- do.call(simfun, SMC_summary$parlist[4,3:7])
save("vergara_sim", file="fitted_sim.rda")

dagan_sim <- do.call(simfun, SMC_summary$parlist[3,3:7])
save("vergara_sim", "dagan_sim", file="fitted_sim.rda")

