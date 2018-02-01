library(emdbook)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/vergara_summ.rda")
load("../data/SMC_vergara.rda")

Nmax <- 10
tolerance <- 0.6

sumlist2 <- simlist2 <- parlist2 <- list()

subyear <- c(1001:1100)
sitesample <- 4

cc <- 1

sigma <- with(as.list(parlist[[4]]),{
    sqrt(2*c(
        var(beta.meanlog),
        var(log(beta.sdlog)),
        var(logit(V)),
        var(logit(epsilon.site)),
        var(log(c_b))
    ))
})

N <- 1
set.seed(101)
while(N <= Nmax) {
    cat(cc, N, "\n")
    pindex <- sample(1:100, 1, prob=ww[[4]])
    pp.sample <- parlist[[4]][pindex,]
    
    pp <- rjump(pp.sample, sigma)
            
    sim <- try(do.call(simfun, c(pp, summarize=FALSE)))
    print(summ <- try(do.call(sumfun, list(sim=sim, subyear=subyear, sitesample=sitesample))))
            
    if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !any(is.na(summ))) {
        print(dist <- sum(abs(vergara_summ - summ)))
        accept <- dist < tolerance
        if (accept) {
            parlist2[N,] <- unlist(pp)
            sumlist2[N,] <- summ
                    
            simlist2[[N]] <- sim
            
            N <- N+1
            save("sumlist2", "simlist2", "parlist2", file="SMC_SHARCNET.rda")
        }
    }
    cc <- cc+1
}

save("sumlist2", "simlist2", "parlist2", file="SMC_SHARCNET.rda")
