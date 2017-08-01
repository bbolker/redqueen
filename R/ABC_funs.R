sumfun <- function(sim, subyear=c(1001:1100),
                   sitesample=4) {
    if (!missing(sitesample)) {
        site <- sample(1:dim(sim$S.count)[2], sitesample)
    } else {
        site <- 1:(dim(S.count)[2])
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
        
        names(summ) <- names(vergara_summ)
        return(summ)
    })
}

simfun <- function(beta.meanlog=1, beta.sdlog=0.5,
                   bU=20, V=0.85,
                   epsilon.site=0.01,
                   n.site=50,
                   subyear=c(1001:1100),
                   summarize=TRUE,
                   discard=TRUE, 
                   sitesample=4, ...) {
    bI <- (1-V)*bU
    beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
    tmax <- max(subyear)
    sim <- three_loci_stochastic_spatial_discrete_model(beta=beta, n.site=n.site, bU=bU, bI=bI, epsilon.site=epsilon.site, tmax=tmax, ...)
    
    if (discard && any(sim$A.count[subyear,] < 0.1)) return(NA)
    
    summary <- sumfun(sim=sim, subyear=subyear, sitesample=sitesample)
    
    if (summarize) {
        return(summary)
    } else {
        return(sim)
    }
}

betafun <- function(beta.meanlog=1, beta.sdlog=0.5,
                    bU=20, V=0.85,
                    epsilon.site=0.01,
                    n.site=4,
                    subyear=c(1001:1100),
                    nsim=10, verbose=FALSE,
                    debug=FALSE, ...) {
    i <- 1
    j <- 1
    combres <- vector('list', nsim)
    while (i <= nsim) {
        res <- simfun(beta.meanlog=beta.meanlog, beta.sdlog=beta.sdlog, 
                      bU=bU, V=V,
                      epsilon.site = epsilon.site,
                      n.site=n.site, subyear=subyear, ...)
        if(!is.na(res[1])) {
            combres[[i]] <- res
            i <- i + 1
        }
        j <- j + 1
        if (debug) {cat(i, ",",  j, "\n"); print(res)}
        if (j > 2*nsim) break
    }
    
    combsum <- colMeans(do.call(rbind, lapply(combres, unlist)))
    if(verbose) print(combsum)
    return(combsum)
}
