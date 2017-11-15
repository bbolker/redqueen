logit <- function(x,eps=1e-10) {
    qlogis(pmax(eps,pmin(1-eps,x)))
}

rprior <- function() {
    list(
        beta.meanlog=rcauchy(1, location=1, scale=2),
        beta.sdlog=rlnorm(1, meanlog=0, sdlog=2),
        V=rbeta(1, shape1=6, shape2=2), ## mean of 0.75
        epsilon.site=rbeta(1, shape1=1, shape2=9), ## mean of 0.1
        n.genotype=rbetabinom(n=1,size=9,prob=3/9,theta=5)+1, ## mean of 3
        c_b=rlnorm(1, meanlog=-0.07, sdlog=0.09)
    )
}

dprior <- function(x) {
    with(as.list(x),{
        dcauchy(beta.meanlog, location=1, scale=2) *
            dlnorm(beta.sdlog, meanlog=0, sdlog=2) *
            dbeta(V, shape1=6, shape2=2) *
            dbeta(epsilon.site, shape1=1, shape2=9) *
            dbetabinom(n.genotype-1, size=9, prob=3/9, theta=5) *
            dlnorm(c_b, meanlog=-0.07, sdlog=0.09)
    })
}

## assuming that all parameters are independent
rjump <- function(x, sigma) {
    with(as.list(x),{
        logit.V <- logit(V)
        logit.epsilon <- logit(epsilon.site)
        
        list(
            beta.meanlog=rnorm(1, mean=beta.meanlog, sd=sigma[1]), 
            beta.sdlog=exp(rnorm(1, mean=log(beta.sdlog), sd=sigma[2])), 
            V=plogis(rnorm(1, mean=logit.V, sd=sigma[3])), 
            epsilon.site=plogis(rnorm(1, mean=logit.epsilon, sd=sigma[4])), 
            n.genotype=rbinom(1, size=9, prob=((n.genotype-0.5)/10))+1, ## avoid p = 0 and 1
            c_b=exp(rnorm(1, mean=log(c_b), sd=sigma[5]))
        )
    })
}

djump <- function(x, theta, sigma) {
    dnorm(x[[1]], mean=theta[[1]], sd=sigma[1]) *
        dnorm(log(x[[2]]), mean=log(theta[[2]]), sd=sigma[2]) *
        dnorm(logit(x[[3]]), logit(theta[[3]]), sigma[3]) *
        dnorm(logit(x[[4]]), mean=logit(theta[[4]]), sigma[4]) *
        dbinom(x[[5]]-1, size=9, prob=(theta[[5]]-0.5)/10) *
        dnorm(log(x[[6]]), mean=log(theta[[6]]), sd=sigma[[5]])
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

simfun <- function(beta.meanlog=1, beta.sdlog=0.5,
                   bU=20, V=0.85,
                   epsilon.site=0.01,
                   n.site=40,
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
