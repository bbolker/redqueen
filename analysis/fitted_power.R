source("../R/powerfun.R")
load("../data/fitted_sim.rda")

sites <- seq(10, 30, by=5)

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
        sub_reslist[[test_name]] <- lapply(sites,
            function(x) powerfun(
                simlist=sim,
                nsite=x,
                nsample=100,
                nsim=10,
                test=test_list[[test_name]]
            )
        )
    }
    simlist[[sim_name]] <- sub_reslist
}
    
save("simlist", file="fitted_power.rda")
