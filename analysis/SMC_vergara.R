library(emdbook)
library(car)
source("../R/util.R")
source("../R/stochastic_model.R")
source("../R/ABC_funs.R")

load("../data/vergara_summ.rda")

rprior <- function() {
    list(
        beta.meanlog=rcauchy(1, location=2, scale=1),
        beta.sdlog=rlnorm(1, meanlog=0, sdlog=1),
        V=rbeta(1, shape1=6, shape2=2), ## mean of 0.75
        epsilon.site=rbeta(1, shape1=1, shape2=99), ## mean of 0.01
        n.genotype=rbetabinom(n=1,size=9,prob=2/9,theta=5)+1, ## mean of 3
        c_b=rlnorm(1, meanlog=-0.1, sdlog=0.1)
    )
}

dprior <- function(x) {
    with(as.list(x),{
        dcauchy(beta.meanlog, location=2, scale=1) *
            dlnorm(beta.sdlog, meanlog=0, sdlog=1) *
            dbeta(V, shape1=6, shape2=2) *
            dbeta(epsilon.site, shape1=1, shape2=99) *
            dbetabinom(n.genotype-1, size=9, prob=2/9, theta=5) *
            dlnorm(c_b, meanlog=-0.1, sdlog=0.1)
    })
}

## assuming that all parameters are independent
rjump <- function(x, sigma) {
    with(as.list(x),{
        logit.V <- car::logit(V, adjust=1e-10)
        logit.epsilon <- car::logit(epsilon.site, adjust=1e-10)
        
        list(
            beta.meanlog=rnorm(1, mean=beta.meanlog, sd=sigma[1]), 
            beta.sdlog=exp(rnorm(1, mean=log(beta.sdlog), sd=sigma[2])), 
            V=plogis(rnorm(1, mean=logit.V, sd=sigma[3])), 
            epsilon.site=plogis(rnorm(1, mean=logit.epsilon, sd=sigma[4])), 
            n.genotype=rbinom(1, size=9, prob=((n.genotype-0.5)/10))+1, ## avoid p = 0 and 1
            c_b=exp(rnorm(1, mean=log(c_b), sd=sigma[5]))
        )
    })
}

djump <- function(x, theta, sigma) {
    dnorm(x[[1]], mean=theta[[1]], sd=sigma[1]) *
        dnorm(log(x[[2]]), mean=log(theta[[2]]), sd=sigma[2]) *
        dnorm(car::logit(x[[3]], adjust=1e-10), mean=car::logit(theta[[3]], adjust=1e-10), sigma[3]) *
        dnorm(car::logit(x[[4]], adjust=1e-10), mean=car::logit(theta[[4]], adjust=1e-10), sigma[4]) *
        dbinom(x[[5]]-1, size=9, prob=(theta[[5]]-0.5)/10) *
        dnorm(log(x[[6]]), mean=log(theta[[6]]), sd=sigma[[5]])
}

Nmax <- 50
tmax <- 3
tolerance <- c(1.8, 0.9, 0.6)

ww <- matrix(NA, ncol=tmax, nrow=Nmax)
sumlist <- simlist <- parlist <- vector("list", tmax)

subyear <- c(1001:1100)
sitesample <- 4

for(t in 1:tmax) {
    parlist[[t]] <- as.data.frame(matrix(NA, ncol=length(rprior()), nrow=Nmax))
    names(parlist[[t]]) <- names(rprior())
    
    sumlist[[t]] <- as.data.frame(matrix(NA, ncol=length(vergara_summ), nrow=Nmax))
    names(sumlist[[t]]) <- names(vergara_summ)
    
    simlist[[t]] <- vector('list', Nmax)

    N <- 1
    if (t== 1) {
        while(N <= Nmax) {
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
                    simlist[[t]][[N]] <- sim
                    ww[N,t] <- 1/Nmax
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
        while(N <= Nmax) {
            cat(t, N, "\n")
            pindex <- sample(1:Nmax, 1, prob=ww[,t-1])
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
                    simlist[[t]][[N]] <- sim
                    
                    ww[N,t] <- dprior(pp)/
                        sum(ww[,t-1] * apply(parlist[[t-1]], 1, djump, theta=pp, sigma=sigma))
                    N <- N+1
                    save("ww", "sumlist", "simlist", "parlist", file="SMC_vergara.rda")
                }
            }
        }
        ww[,t] <- ww[,t]/sum(ww[,t])
    }
}
