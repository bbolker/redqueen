library(dplyr)
source("../R/powerfun.R")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn_sim <- paste0("post_dagan_sim_", batch_num, ".rda")

load(fn_sim)

fn_save <- paste0("lm_power_dagan_", batch_num, ".rda")

site_vec <- seq(10, 40, by=10)
sample_vec <- seq(25, 150, by=25)
nsim <- 100

reslist <- vector('list', 10)

for (i in 1:10) {
    sim <- sim_reslist[[i]]
    
    sitelist <- vector('list', length(site_vec))
    
    for (j in 1:length(site_vec)) {
        nsite <- site_vec[j]
        
        samplelist <- vector('list', length(sample_vec))
        
        for (k in 1:length(sample_vec)) {
            nsample <- sample_vec[k]
            
            print(paste(i, j, k, sep=", "))
            
            target.gen <- 1099:1100
            
            power <- powerfun(sim, 
                              nsim=nsim,
                              nsample=nsample,
                              nsite=nsite,
                              transform="arcsin",
                              target.gen=target.gen,
                              test=list(lm=test_lm))
            
            powerdf <- power
            
            powerdf$sample <- nsample
            
            samplelist[[k]] <- powerdf
        }
        
        sampledf <- do.call("rbind", samplelist)
        
        sampledf$site <- nsite
        
        sitelist[[j]] <- sampledf
    }
    
    sitedf <- do.call("rbind", sitelist)
    
    sitedf$param <- i + batch_num * 10
    
    reslist[[i]] <- sitedf
}

save("reslist", file=fn_save)
