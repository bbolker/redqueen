source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/simfun.R")

L <- load("../data/SMC_summary.rda")

vergara_sim <- do.call(simfun, SMC_summary$parlist[4,3:8])
save("vergara_sim", file="fitted_sim.rda")

dagan_sim <- do.call(simfun, SMC_summary$parlist[3,3:8])
save("vergara_sim", "dagan_sim", file="fitted_sim.rda")

