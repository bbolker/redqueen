library(dplyr)
source("../R/powerfun.R")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn_sim <- paste0("post_dagan_sim_", batch_num, ".rda")

load(fn_sim)

fn_save <- paste0("power_dagan_", batch_num, ".rda")

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
            
            powerlist <- vector('list', 99)
            
            for (gg in 0:98) {
                print(paste(i, j, k, gg, sep=", "))
                
                target.gen <- 1001:1002 + gg
                
                power <- powerfun(sim, 
                         nsim=nsim,
                         nsample=nsample,
                         nsite=nsite,
                         transform="raw",
                         target.gen=target.gen)
                
                power$gen <- target.gen[1]
                
                powerlist[[gg+1]] <- power
            }
            
            powerdf <- do.call("rbind", powerlist)
            
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
