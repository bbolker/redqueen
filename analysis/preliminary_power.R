library(dplyr)
library(tidyr)
load("../data/powersim.rda")

s <- 0.5

target.gen <- 1000

samplesize <- c(20, 50, 100, 200)

for(i in 1:dim(simlist)[2]) {
    for (j in 1:dim(simlist)[1]) {
        sim <- simlist[j,i][[1]]
        sim$S.count[target.gen,]
    }
}


spatial_power <- function(data,
         n.site=18,
         n.sample=100,
         target.gen=2500,
         n.replicate=100) {

    powerres <- vector('list', n.replicate)
    
    for (i in 1:n.replicate) {
        which.site <- sample(1:50, n.site)
        
        spreadres <- data %>% 
            lapply(function(x) data.frame(site=which.site, count=x[target.gen,which.site])) %>%
            bind_rows(.id="type") %>%
            spread(type, count)
        
        sampledf <- spreadres %>% 
            apply(1, function(x) rmultinom(n=1, size=n.sample, prob=x[-1]/sum(x[-1]))) %>%
            t %>%
            as.data.frame %>%
            setNames(names(spreadres)[-1]) %>%
            bind_cols(data.frame(site=which.site))
        
        maledf <- data.frame(
            male=sampledf$MI.count+sampledf$MU.count,
            female=n.sample-(sampledf$MI.count+sampledf$MU.count)
        )
        
        infdf <- data.frame(
            infected=(sampledf$AI.count+sampledf$FI.count+sampledf$MI.count),
            uninfected=n.sample-(sampledf$AI.count+sampledf$FI.count+sampledf$MI.count)
        )
        
        chi_male <- chisq.test(maledf)
        chi_infection <- chisq.test(infdf)
        corr <- cor.test(infdf$infected/n.sample, maledf$male/n.sample, method="spearman")
        
        powerres[[i]] <- data.frame(
            chi_male=chi_male$p.value,
            chi_infection=chi_infection$p.value,
            corr=corr$p.value
        )
        
    }
    
    power <- colSums(do.call("rbind", powerres) < 0.05)/n.replicate
    return(power)
}

sample_vec <- c(20, 50, 100, 500, 1000)
site_vec <- c(10, 20, 30, 40, 50)

(power_list_sample <- lapply(sample_vec, function(x) spatial_power(subres, n.sample=x)))
(power_list_site <- lapply(site_vec, function(x) spatial_power(subres, n.site=x)))
