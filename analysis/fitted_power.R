library(dplyr)
source("../R/powerfun.R")
load("../data/SMC_summary.rda")

sites <- seq(10, 30, by=5)
samples <- seq(25, 150, by=25)
nsim <- 100

simlist <- list(
    dagan=simlist$dagan[[3]],
    vergara=simlist$vergara[[3]],
    mckone=simlist$mckone[[3]]
)

test_list <- list(
    linear=test_lm,
    quadratic=test_quad,
    spearman=test_spearman
)

reslist <- list()

for(sim_name in names(simlist)) {
    sim <- simlist[[sim_name]]
    sub_reslist <- list()
    
    sample_reslist <- vector('list', length=length(samples))
    for(i in 1:length(samples)) {
        res <- lapply(sites,
            function(x) powerfun(
                simlist=sim,
                nsite=x,
                nsample=samples[i],
                nsim=nsim,
                test=test_list,
                verbose=TRUE
            )
        )
        names(res) <- sites
        rres <- res %>%
            bind_rows(.id="sites")
        
        rres$sim <- rep(rep(1:50, each=3*nsim), length(sites))
        
        rres$samples <- samples[i]
        sample_reslist[[i]] <- rres
    }
    
    sample_res <- do.call('rbind', sample_reslist)
    
    reslist[[sim_name]] <- sample_res
}

save("reslist", file="fitted_power.rda")
