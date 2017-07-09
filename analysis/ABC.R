source("../R/util.R")
source("../R/stochastic_model.R")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf", "psex", "pinf", "psex"
))

simfun <- function(meanlog=1, sdlog=0.5,
                   n.site=4,
                   subyear=c(1001:1100), ...) {
    beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
    tmax <- max(subyear)
    sim <- stochastic_spatial_discrete_model(beta=beta, n.site=n.site, tmax=tmax, ...)
    
    if (any(sim$A.count[subyear,] < 0.1)) return(NA)
    
    summary <- with(sim,{
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
    
    return(summary)
}

betafun <- function(meanlog=1, sdlog=0.5,
                    n.site=4,
                    subyear=c(1001:1100),
                    nsim=20, verbose=FALSE,
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
## This parameter seems ok
## betafun(meanlog=2, sdlog=0.2, epsilon.site=0.001, n.site=4, debug=TRUE)

## set up priors...
set.seed(101)
size <- 1000

param.table <- data.frame(
    meanlog=rlnorm(size),
    sdlog=rlnorm(size, meanlog=-2),
    bU=runif(size, min=10, max=30) ## corresponding to virulence 70%-90%
)

param.keep <- rep(NA, size)
res <- vector('list', size)

for (i in 1:size) {
    print(i)
    print(res[[i]] <- do.call(betafun, param.table[i,]))
    print(param.keep[i] <- sqrt(sum((vergara_CV - res[[i]])^2)) < 0.2)
    save("param.table", "res", "param.keep", file="ABC.rda")
}




