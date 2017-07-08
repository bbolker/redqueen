library(readxl)
library(dplyr)
source("../R/util.R")
source("../R/stochastic_model.R")
vergara <- read_excel("../data/vergara2014.xlsx")

cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara_clean <- vergara %>%
    filter(Ploidy != "UNKNOWN")

vergara_p <- vergara_clean %>%
    group_by(Site, Year) %>%
    ## mean prevalence of infection/sex for each site every year
    summarize(pinf=mean(Microphallus), psex=sum(Ploidy=="sexual")/length(Ploidy))

vergara_site_CV <- vergara_p  %>% 
    group_by(Site) %>%
    ## CV across year
    summarize_each(funs(sd(.)/mean(.)), pinf, psex) %>%
    ## mean of CV across site
    select(-Site) %>%
    summarize_all(mean) %>%
    unlist

vergara_year_CV <- vergara_p  %>% 
    group_by(Year) %>%
    ## CV across site
    summarize_each(funs(sd(.)/mean(.)), pinf, psex) %>%
    ## mean of CV across year
    select(-Year) %>%
    summarize_all(mean) %>%
    unlist


simfun <- function(meanlog=1, sdlog=0.5,
                   n.site=4,
                   subyear=c(1001:1100), ...) {
    beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
    tmax <- max(subyear)
    sim <- stochastic_spatial_discrete_model(beta=beta, n.site=n.site, tmax=tmax, ...)
    
    if (any(sim$A.count[1000:tmax,] < 0.1) || any(sim$S.count[1000:tmax,] < 0.1)) return(NA)
    
    summary <- with(sim,{
        N.count <- S.count + A.count
        sl <- list(
            pinf = (SI.count[subyear,] + AI.count[subyear,])/N.count[subyear,],
            psex = S.count[subyear,]/N.count[subyear,]
        )
        
        sl2 <- vector('list', 2)
        
        for (i in 1:2) sl2[[i]] <- lapply(sl, function(x) apply(x, i, function(y) sd(y)/mean(y)))
        
        sl2 <- lapply(sl2, function(x) unlist(lapply(x, mean)))
        names(sl2) <- c("persite", "peryear")
        return(sl2)
    })
    
    return(summary)
}

betafun <- function(meanlog=1, sdlog=0.5,
                    n.site=4,
                    subyear=c(1001:1100),
                    nsim=10, verbose=TRUE,
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
        if(debug) cat(i, ",",  j, "\n"); print(res)
    }
    
    combsum <- colMeans(do.call(rbind, lapply(combres, unlist)))
    if(verbose) print(combsum)
}

betafun(meanlog=2, sdlog=0.1, epsilon.site=0.0001, n.site=4, debug=TRUE)

ss <- stochastic_spatial_discrete_model(beta=rlnorm(4, meanlog=2, sdlog=0.1), n.site=4, epsilon.site=0.0001)
matplot(ss$S.count/(ss$S.count+ss$A.count), type="l")
matplot(ss$A.count, type="l")


