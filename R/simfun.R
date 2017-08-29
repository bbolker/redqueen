sumfun <- function(sim, subyear=c(1001:1100),
                   sitesample=4) {
    if (!missing(sitesample)) {
        site <- sample(1:dim(sim$S.count)[2], sitesample)
    } else {
        site <- 1:(dim(sim$S.count)[2])
    }
    
    with(sim,{
        N.count <- S.count + A.count
        sl <- list(
            pinf = (SI.count[subyear,site] + AI.count[subyear,site])/N.count[subyear,site],
            psex = S.count[subyear,site]/N.count[subyear,site]
        )
        
        sl2 <- vector('list', 2)
        
        for (i in 1:2) sl2[[i]] <- lapply(sl, function(x) apply(x, i, mean))
        
        cv <- unlist(lapply(sl2, function(x) sapply(x, function(y) sd(y)/mean(y))))
        
        mean <- unlist(lapply(sl2, function(x) lapply(x, mean))[[1]])
        
        summ <- c(cv, mean)
        
        names(summ) <- c("pinf.timeCV", "psex.timeCV", "pinf.siteCV", "psex.siteCV", "pinf.mean", "psex.mean")
        return(summ)
    })
}

simfun <- function(nsim=3, 
                   beta.meanlog,
                   beta.sdlog,
                   epsilon.site,
                   n.genotype,
                   V, 
                   c_b,
                   n.site=30) {
    simlist <- vector("list", nsim)
    j <- 1
    while(j <= nsim) {
        beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
        bU <- 20
        bI <- bU * (1-V)
        sim <- stochastic_spatial_discrete_lim_model(
            beta=beta,
            n.site=n.site,
            n.genotype=round(n.genotype),
            epsilon.site=epsilon.site,
            bU=bU,
            bI=bI,
            c_b=c_b
        )
        
        summary <- sumfun(sim=sim, subyear=1000:1100, sitesample=n.site)
        
        if(summary[["psex.mean"]] > 0.001) {
            simlist[[j]] <- sim
            j <- j + 1
        }
        cat(j, "\n")
    }
    simlist
}