library(dplyr)
source("../R/powerfun.R")
load("../data/SMC_summary.rda")

gens <- 1001:1002

sites <- seq(10, 40, by=10)
samples <- seq(25, 150, by=25)
nsim <- 100
ngen <- c(1, 2, 5, 10, 20)

simlist2 <- list(
    vergara=simlist$vergara[[4]],
    mckone=simlist$mckone[[4]]
)

test_list <- list(
    spearman=test_spearman
)

reslist <- vector('list', length(simlist2))
names(reslist) <- names(simlist2)

for(sim_name in names(simlist2)) {
    print(sim_name)
    sim <- simlist2[[sim_name]][1:40] ## TODO: Change this once simulations are over!!!
    
    sample_reslist <- vector('list', length=length(samples))
    for(i in 1:length(samples)) {
        print(i)
        set.seed(i)
        
        res <- lapply(sites,
            function(x) powerfun(
                simlist=sim,
                nsite=x,
                nsample=samples[i],
                nsim=nsim,
                test=test_list,
                transform="raw",
                target.gen=gens,
                verbose=FALSE
            )
        )
        names(res) <- sites
        rres <- res %>%
            bind_rows(.id="sites")
        
        rres$sim <- rep(rep(1:40, each=length(test_list)*nsim), length(sites))
        
        rres$samples <- samples[i]
        sample_reslist[[i]] <- rres
    }
    
    sample_res <- do.call('rbind', sample_reslist)
    
    reslist[[sim_name]] <- sample_res
}

save("reslist", file="fitted_power.rda")
