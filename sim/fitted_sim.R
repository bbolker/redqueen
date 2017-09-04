source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/simfun.R")

L <- load("../data/SMC_summary.rda")

vergara_par <- comb_smc$vergara$parlist[[2]]
dagan_par <- comb_smc$dagan$parlist[[2]]
mckone_par <- comb_smc$mckone$parlist[[2]]

vergara_sim <- dagan_sim <- mckone_sim <- vector('list', nrow(vergara_par))

for(i in 1:nrow(vergara_par)) {
    cat("vergara", i, "\n")
    vergara_sim[[i]] <- do.call(simfun, vergara_par[i,])
    save("vergara_sim", file="fitted_sim.rda")
}

for(i in 1:nrow(dagan_par)) {
    cat("dagan", i, "\n")
    dagan_sim[[i]] <- do.call(simfun, dagan_par[i,])
    save("vergara_sim", "dagan_sim", file="fitted_sim.rda")
}

for(i in 1:nrow(mckone_par)) {
    cat("mckone", i, "\n")
    mckone_sim[[i]] <- do.call(simfun, mckone_par[i,])
    save("vergara_sim", "dagan_sim", "mckone_sim", file="fitted_sim.rda")
}