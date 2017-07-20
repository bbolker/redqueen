source("../R/util.R")
source("../R/stochastic_model.R")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf", "psex", "pinf", "psex"
    ))

sumfun <- function(sim, subyear=c(1001:1100)) {
    with(sim,{
        N.count <- S.count + A.count
        sl <- list(
            pinf = (SI.count[subyear,] + AI.count[subyear,])/N.count[subyear,],
            psex = S.count[subyear,]/N.count[subyear,]
        )
        
        sl2 <- vector('list', 2)
        
        for (i in 1:2) sl2[[i]] <- lapply(sl, function(x) apply(x, i, mean))
        
        sl2 <- lapply(sl2, function(x) unlist(lapply(x, function(y) sd(y)/mean(y))))
        names(sl2) <- c("peryear", "persite")
        return(sl2)
    })
}

simfun <- function(meanlog=1, sdlog=0.5,
                   n.site=4,
                   subyear=c(1001:1100),
                   summarize=TRUE, ...) {
    beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
    tmax <- max(subyear)
    sim <- stochastic_spatial_discrete_model(beta=beta, n.site=n.site, tmax=tmax, ...)
    
    if (any(sim$A.count[subyear,] < 0.1)) return(NA)
    
    summary <- sumfun(sim, subyear)
    
    if (summarize) {
        return(summary)
    } else {
        return(sim)
    }
    
}

betafun <- function(meanlog=1, sdlog=0.5,
                    n.site=4,
                    subyear=c(1001:1100),
                    nsim=10, verbose=FALSE,
                    debug=FALSE,
                    ...) {
    i <- 1
    j <- 1
    combres <- vector('list', nsim)
    while (i <= nsim) {
        res <- simfun(meanlog=meanlog, sdlog=sdlog, n.site=n.site, subyear=subyear, ...)
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

## set up priors...
set.seed(101)
size <- 10000
tol <- 0.4

param.table <- data.frame(
    meanlog=rlnorm(size),
    sdlog=rlnorm(size, meanlog=-2, sdlog=0.5),
    bU=runif(size, min=10, max=30), 
    bI=runif(size, min=1, max=5), ## corresponds to virulence of 50%-97%
    epsilon.site=exp(runif(size, min=log(1e-4), max=log(1e-1)))
)

## exmample
## rr <- do.call(simfun, c(param.table[1000,], summarize=FALSE))
## sumfun(rr)

param.keep <- rep(NA, size)
res <- vector('list', size)

for (i in 1:size) {
    print(i)
    print(res[[i]] <- do.call(betafun, param.table[i,]))
    print(param.keep[i] <- sqrt(sum((vergara_CV - res[[i]])^2)) < tol)
    save("param.table", "res", "param.keep", file="ABC2.rda")
}
