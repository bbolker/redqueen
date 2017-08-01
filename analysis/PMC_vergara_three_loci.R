library(mvtnorm)
source("../R/three_loci_util.R")
source("../R/three_loci_stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/vergara_summ.rda")

rprior <- function() {
    list(
        beta.meanlog=rcauchy(1, location=2, scale=1),
        log.beta.sdlog=rcauchy(1, location=0, scale=0.5),
        logit.V=rcauchy(1, location=1, scale=0.5)
    )
}

dprior <- function(x) {
    with(as.list(x),{
        dcauchy(beta.meanlog, location=2, scale=1) *
            dcauchy(log.beta.sdlog, location=0, scale=0.5) *
            dcauchy(logit.V, location=1, scale=0.5)
    })
}

tpars <- function(par) {
    with(as.list(par),{
        list(beta.meanlog=beta.meanlog,
             beta.sdlog=exp(log.beta.sdlog),
             V=plogis(logit.V))
    })
}

Nmax <- 100
tmax <- 4
tolerance <- c(0.8, 0.4, 0.2, 0.1)

ww <- matrix(NA, ncol=tmax, nrow=Nmax)
sumlist <- tparlist <- parlist <- vector("list", tmax)


for(t in 1:tmax) {
    tparlist[[t]] <- parlist[[t]] <- as.data.frame(matrix(NA, ncol=3, nrow=Nmax))
    names(parlist[[t]]) <- names(rprior())
    names(tparlist[[t]]) <- names(tpars(rprior()))
    
    sumlist[[t]] <- as.data.frame(matrix(NA, ncol=length(vergara_summ), nrow=Nmax))
    names(sumlist[[t]]) <- names(vergara_summ)

    N <- 1
    if (t== 1) {
        while(N <= Nmax) {
            cat(t, N, "\n")
            pp <- rprior()
            tpp <- tpp2 <- tpars(pp)
            
            print(summ <- try(do.call(simfun, tpp2)))
            if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !is.na(summ[1])) {
                print(dist <- sqrt(sum((vergara_summ - summ)^2)))
                accept <- dist < tolerance[t]
                if (accept) {
                    parlist[[t]][N,] <- unlist(pp)
                    tparlist[[t]][N,] <- unlist(tpp)
                    sumlist[[t]][N,] <- summ
                    ww[N,t] <- 1/Nmax
                    N <- N+1
                    save("ww", "sumlist", "tparlist", "parlist", file="PMC_vergara_three_loci.rda")
                }
            }
        }
    } else {
        sigma <- 2*var(parlist[[t-1]])
        while(N <= Nmax) {
            cat(t, N, "\n")
            pindex <- sample(1:Nmax, 1, prob=ww[,t-1])
            mean <- unlist(parlist[[t-1]][pindex,])
            
            pp <- as.list(as.data.frame(rmvnorm(1, mean=mean, sigma=sigma)))
            tpp <- tpp2 <- tpars(pp)
            
            print(summ <- try(do.call(simfun, tpp2)))
            if (!any(is.nan(summ)) && !inherits(summ, "try-error") && !is.na(summ[1])) {
                print(dist <- sqrt(sum((vergara_summ - summ)^2)))
                accept <- dist < tolerance[t]
                if (accept) {
                    parlist[[t]][N,] <- unlist(pp)
                    tparlist[[t]][N,] <- unlist(tpp)
                    sumlist[[t]][N,] <- summ
                    ww[N,t] <- dprior(pp)/sum(ww[,t-1] * dmvnorm(parlist[[t-1]], mean=unlist(pp), sigma=sigma))
                    N <- N+1
                    save("ww", "sumlist", "tparlist", "parlist", file="PMC_vergara_three_loci.rda")
                }
            }
        }
        ww[,t] <- ww[,t]/sum(ww[,t])
    }
}
