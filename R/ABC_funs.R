logit <- function(x,eps=1e-10) {
    qlogis(pmax(eps,pmin(1-eps,x)))
}

rprior <- function() {
    list(
        beta.mean=rgamma(1, shape=2, scale=10),
        beta.cv=rgamma(1, shape=2, scale=0.5),
        V=rbeta(1, shape1=6, shape2=2), ## mean of 0.75
        epsilon.site=rbeta(1, shape1=1, shape2=9), ## mean of 0.1
        n.genotype=rbetabinom(n=1,size=9,prob=3/9,theta=5)+1, ## mean of 3
        c_b=rlnorm(1, meanlog=-0.07, sdlog=0.09)
    )
}

dprior <- function(x) {
    with(as.list(x),{
        dgamma(beta.mean, shape=2, scale=10) *
            dgamma(beta.cv, shape=2, scale=0.5) *
            dbeta(V, shape1=6, shape2=2) *
            dbeta(epsilon.site, shape1=1, shape2=9) *
            dbetabinom(round(n.genotype)-1, size=9, prob=3/9, theta=5) *
            dlnorm(c_b, meanlog=-0.07, sdlog=0.09)
    })
}

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

simfun <- function(beta.mean=10, beta.cv=0.5,
                   bU=20, V=0.85,
                   epsilon.site=0.01,
                   n.site=40,
                   n.genotype=1,
                   c_b=1,
                   subyear=c(1001:1100),
                   summarize=TRUE,
                   discard=TRUE, 
                   sitesample=4, ...) {
    n.genotype <- round(n.genotype)
    
    bI <- (1-V)*bU
    
    shape <- 1/beta.cv^2
    scale <- beta.mean/shape
    
    beta <- rgamma(n.site, shape=shape, scale=scale)
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
