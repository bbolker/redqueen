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

simfun <- function(beta.meanlog=1, beta.sdlog=0.5,
                   bU=20, V=0.85,
                   epsilon.site=0.01,
                   n.site=30,
                   n.genotype=1,
                   c_b=1,
                   subyear=c(1001:1100),
                   summarize=TRUE,
                   discard=TRUE, 
                   sitesample=4, ...) {
    bI <- (1-V)*bU
    beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
    tmax <- max(subyear)
    sim <- stochastic_spatial_discrete_lim_model(beta=beta, 
                                                 n.site=n.site, 
                                                 bU=bU, bI=bI, 
                                                 epsilon.site=epsilon.site,
                                                 c_b=c_b,
                                                 n.genotype=n.genotype,
                                                 tmax=tmax, ...)
    
    summary <- sumfun(sim=sim, subyear=subyear, sitesample=sitesample)
    
    if(discard && summary[["psex.mean"]] < 0.001)
        return(NA)
    
    if (summarize) {
        return(summary)
    } else {
        return(sim)
    }
}
