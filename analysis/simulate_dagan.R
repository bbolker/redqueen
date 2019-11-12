library(emdbook)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")
source("../R/powerfun.R")

load("../data/SMC_dagan.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])

fn_sim <- paste0("post_dagan_sim_", batch_num, ".rda")
fn_summ <- paste0("post_dagan_summ_", batch_num, ".rda")
fn_data <- paste0("post_dagan_data_", batch_num, ".rda")

param <- parlist[[4]][1:10+10*batch_num,]

gen <- subyear <- c(1001:1100)
sitesample <- 22

n.site <- 40
bU <- 20

nsim <- 10
nsum <- 100

sim_reslist <- summ_reslist  <- simdata_reslist <- vector('list', nrow(param))

set.seed(101)
for (i in 1:nrow(param)) {
    
    sim_subreslist <- summ_subreslist <- simdata_subreslist <- vector('list', nrow(param))
    
    j <- 1
    
    while (j <= nsim) {
        print(paste(i, j, sep=","))
        
        sim <- try(do.call(simfun, c(param[i,], summarize=FALSE, discard=FALSE)))
        
        if (!inherits(sim, "try-error") && !is.na(sim)) {
            sumdf <- as.data.frame(t(replicate(nsum, sumfun(sim, sitesample=sitesample))))
            sumdf$weight <- ww[[3]][i+10*batch_num]
            sumdf$sim <- j
            
            simdata <- with(sim,{
                N.count <- S.count + A.count
                S <- (S.count/N.count)[gen,]
                I <- ((SI.count + AI.count)/N.count)[gen,]
                data.frame(
                    infected=colMeans(I),
                    sexual=colMeans(S)
                )
            })
            
            sumdf$sim <- j
            
            sim_subreslist[[j]] <- sim
            summ_subreslist[[j]] <- sumdf
            simdata_subreslist[[j]] <- simdata
            
            j <- j + 1
        }
    }
    
    sim_reslist[[i]] <- sim_subreslist
    
    tempdf <- do.call("rbind", summ_subreslist)
    tempdf$param <- i+10*batch_num
    
    summ_reslist[[i]] <- tempdf
    
    tempdf2 <- do.call("rbind", simdata_subreslist)
    tempdf2$param <- i+10*batch_num
    
    simdata_reslist[[i]] <- tempdf2
    
    save("sim_reslist", file=fn_sim)
    save("summ_reslist", file=fn_summ)
    save("simdata_reslist", file=fn_data)
}

save("sim_reslist", file=fn_sim)
save("summ_reslist", file=fn_summ)
save("simdata_reslist", file=fn_data)
