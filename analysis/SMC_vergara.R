library(emdbook)
library(car)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/vergara_summ.rda")

Nmax <- c(25, 50, 100)
tmax <- 3
tolerance <- c(1.8, 0.9, 0.6)

ww <- vector('list', 3)
sumlist <- simlist <- parlist <- vector("list", tmax)

subyear <- c(1001:1100)
sitesample <- 4

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
            cat(t, N, "\n")
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
        }
    } else {
        sigma <- with(as.list(parlist[[t-1]]),{
            sqrt(2*c(
                var(beta.meanlog),
                var(log(beta.sdlog)),
                var(car::logit(V, adjust=1e-10)),
                var(car::logit(epsilon.site, adjust=1e-10)),
                var(log(c_b))
            ))
        })
        while(N <= Nmax[t]) {
            cat(t, N, "\n")
            pindex <- sample(1:Nmax[t-1], 1, prob=ww[[t-1]])
            pp.sample <- parlist[[t-1]][pindex,]
            
            pp <- rjump(pp.sample, sigma)
            
            sim <- try(do.call(simfun, c(pp, summarize=FALSE)))
            print(summ <- try(do.call(sumfun, list(sim=sim, subyear=subyear, sitesample=sitesample))))
            
            if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !any(is.na(summ))) {
                print(dist <- sum(abs(vergara_summ - summ)))
                accept <- dist < tolerance[t]
                if (accept) {
                    parlist[[t]][N,] <- unlist(pp)
                    sumlist[[t]][N,] <- summ
                    
                    if (t == tmax) simlist[[t]][[N]] <- sim
                    
                    ww[[t]][N] <- dprior(pp)/
                        sum(ww[[t-1]] * apply(parlist[[t-1]], 1, djump, theta=pp, sigma=sigma))
                    
                    N <- N+1
                    save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")
                }
            }
        }
        ww[[t]] <- ww[[t]]/sum(ww[[t]])
    }
}

save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")