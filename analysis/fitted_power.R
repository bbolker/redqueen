source("../R/powerfun.R")
load("../data/fitted_sim.rda")

sites <- seq(10, 30, by=5)
samples <- seq(25, 200, by=25)
nsim <- 10

simlist <- list(
    dagan=dagan_sim,
    vergara=vergara_sim
)

test_list <- list(
    linear=test_lm,
    quadratic=test_quad,
    spearman=test_spearman
)

reslist <- vector('list', length(simlist))

for(sim_name in names(simlist)) {
    cat(sim_name)
    sim <- simlist[[sim_name]]
    sub_reslist <- list()
    for(test_name in names(test_list)) {
        cat(test_name)
        sample_reslist <- vector('list', length=length(samples))
        for(i in 1:length(samples)) {
            res <- lapply(sites,
                function(x) powerfun(
                    simlist=sim,
                    nsite=x,
                    nsample=samples[i],
                    nsim=nsim,
                    test=test_list[[test_name]]
                )
            )
            bres <- do.call('rbind', res)
            bres$sites <- rep(sites, each=nsim*length(sim))
                
            sample_reslist[[i]] <- bres
        }
        
        sample_res <- do.call('rbind', sample_reslist)
        sample_res$samples <- rep(samples, each=nsim*length(sim)*length(sites))
        
        sub_reslist[[test_name]] <- sample_res
    }
    simlist[[sim_name]] <- sub_reslist
}
    
save("simlist", file="fitted_power.rda")
