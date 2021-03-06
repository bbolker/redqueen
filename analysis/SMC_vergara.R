library(emdbook)
library(MASS)
library(Matrix)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/vergara_summ.rda")

Nmax <- c(100, 100, 100, 100)
tmax <- 4
tolerance <- c(2.4, 1.2, 0.9, 0.6)

ww <- vector('list', tmax )
sumlist <- simlist <- parlist <- vector("list", tmax)

subyear <- c(1001:1100)
sitesample <- 4

cc <- 1

set.seed(101)
for(t in 1:tmax) {
    ww[[t]] <- rep(NA, Nmax[t])
    
    parlist[[t]] <- as.data.frame(matrix(NA, ncol=length(rprior()), nrow=Nmax[t]))
    names(parlist[[t]]) <- names(rprior())
    
    sumlist[[t]] <- as.data.frame(matrix(NA, ncol=length(vergara_summ), nrow=Nmax[t]))
    names(sumlist[[t]]) <- names(vergara_summ)
    
    simlist[[t]] <- vector('list', Nmax[t])

    N <- 1
    if (t== 1) {
        while(N <= Nmax[t]) {
            cat(t, cc, N, "\n")
            pp <- rprior()
            
            sim <- try(do.call(simfun, c(pp, summarize=FALSE)))
            print(summ <- try(do.call(sumfun, list(sim=sim, subyear=subyear, sitesample=sitesample))))
            
            if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !any(is.na(summ))) {
                print(dist <- sum(abs(vergara_summ - summ)))
                accept <- dist < tolerance[t]
                if (accept) {
                    parlist[[t]][N,] <- unlist(pp)
                    sumlist[[t]][N,] <- summ
                    ww[[t]][N] <- 1/Nmax[t]
                    N <- N+1
                    save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")
                }
            }
            cc <- cc+1
        }
    } else {
        Sigma <- as.matrix(nearPD(2 * cov.wt(parlist[[t-1]], wt=ww[[t-1]])$cov)[[1]])
        
        while(N <= Nmax[t]) {
            cat(t, cc, N, "\n")
            
            get.sample <- FALSE
            
            while (!get.sample) {
                pindex <- sample(1:Nmax[t-1], 1, prob=ww[[t-1]])
                pp.sample <- parlist[[t-1]][pindex,]
                
                pp <-  mvrnorm(1, mu=unlist(pp.sample), Sigma=Sigma)
                
                suppressWarnings(get.sample <- (dprior(pp) > 0) && !is.nan(dprior(pp)))
            }
            
            sim <- try(do.call(simfun, c(as.list(pp), summarize=FALSE)))
            print(summ <- try(do.call(sumfun, list(sim=sim, subyear=subyear, sitesample=sitesample))))
            
            if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !any(is.na(summ))) {
                print(dist <- sum(abs(vergara_summ - summ)))
                accept <- dist < tolerance[t]
                if (accept) {
                    parlist[[t]][N,] <- pp
                    sumlist[[t]][N,] <- summ
                    
                    if (t == tmax) simlist[[t]][[N]] <- sim
                    
                    ww[[t]][N] <- dprior(pp)/
                        sum(ww[[t-1]] * apply(parlist[[t-1]], 1, dmvnorm, mu=pp, Sigma=Sigma))
                    
                    N <- N+1
                    save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")
                }
            }
            cc <- cc+1
        }
        ww[[t]] <- ww[[t]]/sum(ww[[t]])
    }
}

save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")