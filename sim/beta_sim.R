source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/simfun.R")

L <- load("../data/SMC_summary.rda")

verg_par <- SMC_summary$parlist[4,-c(1:2)]

nsim <- 20

meanlog <- seq(-0.5, 3.5, by=1)

simlist <- vector("list", length(meanlog))
names(simlist) <- meanlog

for(i in 1:length(meanlog)) {
    cat(meanlog[[i]], "\n")
    pp <- verg_par
    pp$beta.meanlog <- meanlog[i]
    simlist[[i]] <- do.call(simfun, c(pp, nsim=nsim))
    save("simlist", file="beta_sim.rda")
}
