library(emdbook)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/powerfun.R")

load("../data/SMC_vergara.rda")

meanpar <- colMeans(parlist[[4]])

n.site <- 40
bU <- 20
bI <- (1-meanpar[["V"]])*bU
epsilon.site <- meanpar[["epsilon.site"]]
c_b <- meanpar[["c_b"]]
n.genotype <- round(meanpar[["n.genotype"]])

level <- 0.05

nm <- 8
ns <- 10

beta.meanlog <- seq(0, 2.5, length=nm)
beta.sdlog <- seq(0.2, 2, length=ns)

q <- (2*(1:n.site)-1)/(2*n.site)

reslist <- vector('list', nm * ns)

set.seed(101)
for (i in 1:nm) {
    for (j in 1:ns) {
        print(c(i, j))
        beta <- qlnorm(q, meanlog=beta.meanlog[i], sdlog=beta.sdlog[j])
        sim <- replicate(5, stochastic_spatial_discrete_lim_model(beta=beta, 
                                                     n.site=n.site, 
                                                     bU=bU, bI=bI, 
                                                     epsilon.site=epsilon.site,
                                                     c_b=c_b,
                                                     n.genotype=n.genotype,
                                                     tmax=1100), simplify=FALSE)
        power <- powerfun(sim, nsim=1000, nsite=20, transform='raw', nsample=100, target.gen=1001:1002)
        pp <- sum(power$p.value < level & power$effect.size > 0, na.rm=TRUE)/5000
        print(pp)
        
        reslist[[10*(i-1)+j]] <- data.frame(
            beta.meanlog=beta.meanlog[i],
            beta.sdlog=beta.sdlog[j],
            median.effect=median(power$effect.size, na.rm=TRUE),
            power=pp
        )
        save("reslist", file="param_power.rda")
    }
}

save("reslist", file="param_power.rda")
