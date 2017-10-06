library(dplyr)
source("../R/powerfun.R")
load("../data/SMC_summary.rda")

gens <- 1001

sites <- seq(10, 30, by=10)
samples <- seq(25, 150, by=25)
nsim <- 100

simlist <- list(
    dagan=simlist$dagan[[3]],
    vergara=simlist$vergara[[3]],
    mckone=simlist$mckone[[3]]
)

test_list <- list(
    spearman=test_spearman,
    quadratic_rq=test_quad_rq
)

reslist <- list()

for(sim_name in names(simlist)) {
    print(sim_name)
    sim <- simlist[[sim_name]]
    sub_reslist <- list()
    
    sample_reslist <- vector('list', length=length(samples))
    for(i in 1:length(samples)) {
        print(i)
        set.seed(i)
        gg <- 1001
        
        res <- lapply(sites,
            function(x) powerfun(
                simlist=sim,
                nsite=x,
                nsample=samples[i],
                nsim=nsim,
                test=test_list,
                transform="raw",
                target.gen=gg,
                verbose=TRUE
            )
        )
        names(res) <- sites
        rres <- res %>%
            bind_rows(.id="sites")
        
        rres$sim <- rep(rep(1:50, each=length(test_list)*nsim), length(sites))
        
        rres$samples <- samples[i]
        sample_reslist[[i]] <- rres
    }
    
    sample_res <- do.call('rbind', sample_reslist)
    
    reslist[[sim_name]] <- sample_res
}

save("reslist", file="fitted_power.rda")

