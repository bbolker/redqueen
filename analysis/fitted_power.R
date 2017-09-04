source("../R/powerfun.R")
load("../data/fitted_sim.rda")

sites <- seq(10, 30, by=5)
samples <- seq(25, 150, by=25)
nsim <- 10

simlist <- list(
    dagan=dagan_sim,
    vergara=vergara_sim,
    mckone=mckone_sim
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
        sim_reslist <- vector('list', length=length(sim))
        for(j in 1:length(sim)) {
            res <- try(lapply(sites,
                function(x) powerfun(
                    simlist=sim[[j]],
                    nsite=x,
                    nsample=samples[i],
                    nsim=nsim,
                    test=test_list
                )
            ))
            
            sim_reslist[[j]] <- try(do.call('rbind', res))
            if(!inherits(sim_reslist[[j]], "try-error")) {
                sim_reslist[[j]]$sim <- j
                sim_reslist[[j]]$sites <- rep(sites, each=nsim*length(sim[[j]]))
            } 
            cat(sim_name, i, j, "\n")
        }
            
        sample_reslist[[i]] <- do.call('rbind', sim_reslist)
        sample_reslist[[i]]$samples <- samples[i]
    }
    
    sample_res <- do.call('rbind', sample_reslist)
    
    reslist[[sim_name]] <- sample_res
}

save("reslist", file="fitted_power.rda")
