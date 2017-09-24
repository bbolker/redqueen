source("../R/util.R")
source("../R/stochastic_model.R")

load("../data/SMC_summary.rda")

parlist <- list(
    dagan=SMC_summary$parlist[7,],
    mckone=SMC_summary$parlist[8,],
    vergara=SMC_summary$parlist[9,]
)

beta_vec <- exp(seq(log(1), log(100), by=0.1))
cost_vec <- c(1, 1.1, 1.5, 2, 2.5, 3)

bU <- 20

nsim <- 50

reslist <- vector('list', length(parlist))

target.gen <- 1001:1100

for(j in 1:length(parlist)) {
    rreslist <- vector('list', length(beta_vec))
    
    V <- parlist[[j]]$V
    
    for(beta in beta_vec) {
        subreslist <- vector('list', length(cost_vec))
        for(cost in cost_vec) {
            print(paste(cost, beta, sep=", "))
            c_b <- 2/cost
            sumlist <- vector('list', nsim)
            for(i in 1:nsim) {
                sim <- stochastic_spatial_discrete_lim_model(
                    n.site=1,
                    n.genotype=round(parlist[[j]]$n.genotype),
                    epsilon.site=0,
                    beta=beta,
                    c_b=c_b,
                    bU=20,
                    bI=(1-V)*bU
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
            
            subreslist[[which(cost==cost_vec)]] <- as.data.frame(t(apply(sumdf, 2, mean, na.rm=TRUE)))
        }
        ss <- do.call('rbind', subreslist)
        ss$cost <- cost_vec
        rreslist[[which(beta==beta_vec)]] <- ss
    }
    res <- do.call('rbind', rreslist)
    
    res$beta <- rep(beta_vec, each=length(cost_vec))
    
    reslist[[j]] <- res
}

names(reslist) <- names(parlist)

save("reslist", file="fitted_beta_sim.rda")


