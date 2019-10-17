library(emdbook)
library(MASS)
library(Matrix)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/mckone_summ.rda")
load("../data/SMC_mckone.rda")

argvals <- commandArgs(trailingOnly=TRUE)
batch_num <- as.numeric(argvals[1])
fn <- paste0("SMC_mckone_final_", batch_num, ".rda")

tolerance <- 0.4

mckone_var <- names(mckone_summ)

subyear <- c(1001:1100)
sitesample <- 18

cc <- 1

t <- 4

Sigma <- as.matrix(nearPD(2 * cov.wt(parlist[[t-1]], wt=ww[[t-1]])$cov)[[1]])

N <- 0
Nmax <- 1

set.seed(batch_num)
while(N < Nmax) {
    cat(t, cc, N, "\n")
    
    get.sample <- FALSE
    
    while (!get.sample) {
        pindex <- sample(1:100, 1, prob=ww[[t-1]])
        pp.sample <- parlist[[t-1]][pindex,]
        
        pp <-  mvrnorm(1, mu=unlist(pp.sample), Sigma=Sigma)
        
        suppressWarnings(get.sample <- (dprior(pp) > 0) && !is.nan(dprior(pp)))
    }
    
    sim <- try(do.call(simfun, c(as.list(pp), summarize=FALSE)))
    print(summ <- try(do.call(sumfun, list(sim=sim, subyear=subyear, sitesample=sitesample))))
    
    if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !any(is.na(summ))) {
        print(dist <- sum(abs(mckone_summ - summ[mckone_var])))
        accept <- dist < tolerance
        if (accept) {
            accept_param <- pp
            accept_summ <- summ[mckone_var]
            
            accept_sim <- sim
            
            accept_weight <- dprior(pp)/
                sum(ww[[t-1]] * apply(parlist[[t-1]], 1, dmvnorm, mu=pp, Sigma=Sigma))
            
            N <- N+1
        }
    }
    cc <- cc+1
}

save("accept_param", "accept_summ", "accept_sim", "accept_weight", file=fn)
