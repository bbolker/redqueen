source("../R/three_loci_util.R")
source("../R/three_loci_stochastic_model.R")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf.timeCV", "psex.timeCV", "pinf.siteCV", "psex.siteCV"
    ))

vergara_mean <- structure(
    c(0.433838731993719, 0.703887586155451), 
    .Names = c("pinf.mean", 
               "psex.mean"))

vergara_summ <- c(vergara_CV, vergara_mean)

sumfun <- function(sim, subyear=c(1001:1100)) {
    with(sim,{
        N.count <- S.count + A.count
        sl <- list(
            pinf = (SI.count[subyear,] + AI.count[subyear,])/N.count[subyear,],
            psex = S.count[subyear,]/N.count[subyear,]
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
                   n.site=4,
                   subyear=c(1001:1100),
                   summarize=TRUE, ...) {
    bI <- (1-V)*bU
    beta <- rlnorm(n.site, meanlog=beta.meanlog, sdlog=beta.sdlog)
    tmax <- max(subyear)
    sim <- three_loci_stochastic_spatial_discrete_model(beta=beta, n.site=n.site, bU=bU, bI=bI, epsilon.site=epsilon.site, tmax=tmax, ...)
    
    if (any(sim$A.count[subyear,] < 0.1)) return(NA)
    
    summary <- sumfun(sim, subyear)
    
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

## set up priors...
set.seed(101)
size <- 10000
tol <- 0.4

param.table <- data.frame(
    beta.meanlog=rnorm(size, mean=2, sd=5),
    beta.sdlog=rlnorm(size, meanlog=-2),
    bU=rlnorm(size, meanlog=3, sdlog=0.5), 
    V=rbeta(size, shape1=17, shape2=3), ## virulence
    epsilon.site=exp(runif(size, min=log(1e-4), max=log(1e-1)))
)

## exmample
## rr <- do.call(simfun, c(param.table[1000,], summarize=FALSE))
## sumfun(rr)

param.keep <- rep(NA, size)
res <- vector('list', size)

for (i in 1:size) {
    print(i)
    print(res[[i]] <- try(do.call(betafun, param.table[i,])))
    if(!inherits(res[[i]], "try-error")){
        print(param.keep[i] <- sqrt(sum((vergara_summ - res[[i]])^2)) < tol)
    }
    cat("accepted:", sum(param.keep, na.rm=TRUE), "\n")
    save("param.table", "res", "param.keep", file="ABC_vergara_three_loci.rda")
}
