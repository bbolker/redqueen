source("../R/util.R")
source("../R/stochastic_model.R")

load("../data/SMC_summary.rda")

parlist <- list(
    dagan=SMC_summary$parlist[7,],
    mckone=SMC_summary$parlist[8,],
    vergara=SMC_summary$parlist[9,]
)

bU <- 20

nsite <- nq <- 30

reslist <- vector('list', length(parlist))

target.gen <- 1001:1100

for(j in 1:length(parlist)) {
    V <- parlist[[j]]$V
    q <- (2*(1:nq)-1)/(2*nq)
    beta <- qlnorm(q, meanlog=parlist[[j]]$beta.meanlog, sdlog=parlist[[j]]$beta.sdlog)
    
    sim <- stochastic_spatial_discrete_lim_model(
        n.site=nsite,
        n.genotype=round(parlist[[j]]$n.genotype),
        epsilon.site=parlist[[j]]$epsilon.site,
        beta=beta,
        c_b=parlist[[j]]$c_b,
        bU=bU,
        bI=(1-V)*bU
    )
    
    df <- with(sim,{
        N.count <- S.count + A.count
        S <- (S.count/N.count)[target.gen,]
        I <- ((SI.count + AI.count)/N.count)[target.gen,]
        data.frame(
            infected=colMeans(I),
            sexual=colMeans(S)
        )
    })
    
    reslist[[j]] <- df
}

names(reslist) <- names(parlist)

save("reslist", file="fitted_beta_sim.rda")
