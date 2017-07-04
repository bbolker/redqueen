library(readxl)
library(dplyr)
source("../R/util.R")
source("../R/stochastic_model.R")
vergara <- read_excel("../data/vergara2014.xlsx")

cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara_clean <- vergara %>%
    filter(Ploidy != "UNKNOWN")

vergara_CV <- vergara_clean %>%
    group_by(Site, Year) %>%
    ## mean prevalence of infection/sex for each site every year
    summarize(pinf=mean(Microphallus), psex=sum(Ploidy=="sexual")/length(Ploidy)) %>% 
    ## CV across year
    summarize_each(funs(sd(.)/mean(.)), pinf, psex) %>%
    ## mean of CV across site
    select(-Site) %>%
    summarize_all(mean) %>%
    unlist

simfun <- function(meanlog=1, sdlog=0.5,
                   n.site=4,
                   subyear=c(2001:2100)) {
    beta <- rlnorm(n.site, meanlog=meanlog, sdlog=sdlog)
    tmax <- max(subyear)
    sim <- stochastic_spatial_discrete_model(beta=beta, n.site=n.site, tmax=tmax, bI=0.1)
    
    summary <- with(sim,{
        N.count <- S.count + A.count
        sl <- list(
            pinf = (SI.count[subyear,] + AI.count[subyear,])/N.count[subyear,],
            psex = S.count[subyear,]/N.count[subyear,]
        )
        
        cv_persite <- lapply(sl, function(x) apply(x, 2, function(y) sd(y)/mean(y)))
        unlist(lapply(cv_persite, mean))
    })
    
    return(summary)
}

betafun <- function(meanlog=1, sdlog=0.5,
                    n.site=4,
                    subyear=c(2001:2100),
                    nsim=10, verbose=TRUE) {
    combres <- replicate(nsim, simfun(meanlog=meanlog, sdlog=sdlog, n.site=n.site, subyear=subyear))
    combsum <- rowMeans(combres)
    if(verbos) print(combsum)
}

betafun(meanlog=1.3, sdlog=0.3)

