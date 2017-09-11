source("../R/util.R")
source("../R/stochastic_model.R")

beta_vec <- exp(seq(log(1), log(50), by=0.1))
c_b_vec <- seq(0.5, 1.5, by=0.1)

reslist <- vector('list', length(beta_vec))

target.gen <- 1001:1100

for(beta in beta_vec) {
    subreslist <- vector('list', length(c_b_vec))
    for(c_b in c_b_vec) {
        print(paste(c_b, beta, sep=", "))
        sim <- stochastic_spatial_discrete_lim_model(
            n.site=1,
            n.genotype=1,
            epsilon.site=0,
            beta=beta,
            c_b=c_b,
            seed=101
        )
        subreslist[[which(c_b==c_b_vec)]] <- with(sim,{
            N.count <- S.count+A.count
            I.count <- SI.count+AI.count
            data.frame(
                psex=mean((S.count/N.count)[target.gen]),
                pinf=mean((I.count/N.count)[target.gen])
            )
        })
        
    }
    ss <- do.call('rbind', subreslist)
    ss$c_b <- c_b_vec
    reslist[[which(beta==beta_vec)]] <- ss
}

res <- do.call('rbind', reslist)

res$beta <- rep(beta_vec, each=length(c_b_vec))

save("res", file="beta_sim.rda")


