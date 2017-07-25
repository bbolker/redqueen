source("../R/powerfun.R")
load("../data/vergara_fitted.rda")

## power by sample vs sites

nsamples <- seq(10, 100, by=10)
nsites <- seq(10, 50, by=10)

power_list <- effect_list <- vector('list', length(nsites))
    
for(i in 1:length(nsites)) {
    cat(i)
    res <- lapply(nsamples, function(s) powerfun(simlist, nsample=s, nsite = nsites[[i]]))
    power_list[[i]] <- data.frame(
        nsamples,
        power=sapply(res, "[[", index="power")
    )
    effect_list[[i]] <- cbind(nsamples, do.call("rbind" ,lapply(res, "[[", index="effect")))
    rownames(effect_list[[i]]) <- NULL
    save("effect_list", "power_list", file="vergara_fitted_power.rda")
}

## temporal variation in effect?

gens <- 1001:1100

gen_power_list <- lapply(gens, function(g) powerfun(simlist, nsample=40, nsite = 25, target.gen=g))

gen_power <- data.frame(
    gens,
    power=sapply(gen_power_list, "[[", index="power")
)

gen_effect <- cbind(gens, do.call("rbind" ,lapply(gen_power_list, "[[", index="effect")))
rownames(gen_effect) <- NULL

save("effect_list", "power_list", "gen_power", "gen_effect", file="vergara_fitted_power.rda")

## summary statistics?

vergara_fitted_summary <- lapply(simlist, function(sim) {
    S <- sim$S.count[gens,]
    SI <- sim$SI.count[gens,]
    A <- sim$A.count[gens,]
    AI <- sim$AI.count[gens,]
    N <- S+A
    data.frame(
        psex=rowMeans(S/N),
        pinf=rowMeans((SI+AI)/N),
        psex_inf=rowMeans(SI/S),
        pasex_inf=rowMeans(AI/A)
    )
})

save("vergara_fitted_summary", "effect_list", "power_list", "gen_power", "gen_effect", file="vergara_fitted_power.rda")
