source("../R/util.R")
source("../R/stochastic_model.R")

beta_vec <- exp(seq(log(1), log(50), by=0.1))
cost_vec <- seq(1, 3, by=0.1)

nsim <- 50

reslist <- vector('list', length(beta_vec))

target.gen <- 1001:1100

for(beta in beta_vec) {
    subreslist <- vector('list', length(cost_vec))
    for(cost in cost_vec) {
        print(paste(cost, beta, sep=", "))
        c_b <- 2/cost
        sumlist <- vector('list', nsim)
        for(i in 1:nsim) {
            sim <- stochastic_spatial_discrete_lim_model(
                n.site=1,
                n.genotype=1,
                epsilon.site=0,
                beta=beta,
                c_b=c_b
            )
            sumlist[[i]] <- with(sim,{
                N.count <- S.count+A.count
                I.count <- SI.count+AI.count
                data.frame(
                    psex=mean((S.count/N.count)[target.gen]),
                    pinf=mean((I.count/N.count)[target.gen])
                )
            })
        }
        
        sumdf <- do.call('rbind', sumlist)
        
        subreslist[[which(cost==cost_vec)]] <- as.data.frame(t(apply(sumdf, 2, mean)))
    }
    ss <- do.call('rbind', subreslist)
    ss$cost <- cost_vec
    reslist[[which(beta==beta_vec)]] <- ss
}

res <- do.call('rbind', reslist)

res$beta <- rep(beta_vec, each=length(cost_vec))

save("res", file="beta_sim.rda")


