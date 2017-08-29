source("../R/powerfun.R")
load("../data/fitted_sim.rda")

sites <- seq(10, 30, by=5)
samples <- seq(25, 200, by=25)
nsim <- 10

simlist <- list(
##    dagan=dagan_sim,
    vergara=vergara_sim
)

test_list <- list(
    linear=test_lm,
    quadratic=test_quad,
    spearman=test_spearman
)

reslist <- vector('list', length(simlist))

for(sim_name in names(simlist)) {
    sim <- simlist[[sim_name]]
    sub_reslist <- list()
    for(test_name in names(test_list)) {
        sample_reslist <- vector('list', length=length(samples))
        for(i in 1:length(samples)) {
            sim_reslist <- vector('list', length=length(sim))
            for(j in 1:length(sim)) {
                res <- try(lapply(sites,
                    function(x) powerfun(
                        simlist=sim[[j]],
                        nsite=x,
                        nsample=samples[i],
                        nsim=nsim,
                        test=test_list[[test_name]]
                    )
                ))
                sim_reslist[[j]] <- try(do.call('rbind', res))
                if(!inherits(sim_reslist[[j]], "try-error")) sim_reslist[[j]]$sim <- j
                cat(sim_name, test_name, i, j, "\n")
            }
                
            sample_reslist[[i]] <- do.call('rbind', sim_reslist)
        }
        
        sample_res <- do.call('rbind', sample_reslist)
        sample_res$samples <- rep(samples, each=nsim*length(sim)*length(sites))
        
        sub_reslist[[test_name]] <- sample_res
    }
    simlist[[sim_name]] <- sub_reslist
}
    
save("simlist", file="fitted_power.rda")
