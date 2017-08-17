test_lm <- function(sample) {
    fit <- lm(sexual~infected, data=sample)
    data.frame(
        effect.size=fit$coefficients[2],
        p.value=summary(fit)$coefficients[2,'Pr(>|t|)']
    )
}

test_quad <- function(sample) {
    fit <- lm(sexual~I(infected)+I(infected^2), data=sample)
    data.frame(
        effect.size=fit$coefficients[3],
        p.value=summary(fit)$coefficients[3,'Pr(>|t|)']
    )
}

test_spearman <- function(sample) {
    fit <- cor.test(sample$infected, sample$sexual, method="spearman")
    data.frame(
        effect.size=fit$estimate,
        p.value=fit$p.value
    )
}

sample_sim <- function(sim,
                  nsample,
                  nsite,
                  transform=c("logit", "arcsin", "raw"),
                  target.gen) {
    transform <- match.arg(transform)
    transfun <- switch(transform,
        logit={
            function(p) {
                if(any(p==0) || any(p==1)) {
                    log((p+1/nsample)/(1-p+1/nsample))
                } else {
                    log(p/(1-p))
                }
            } 
        },
        arcsin={
            function(p) arcsin(sqrt(p))
        },
        raw={
            function(p) p
        }
    )
    
    sites <- ncol(sim$S.count)
    target.sites <- sample(1:sites, nsite)
    
    S <- sim$S.count[target.gen, target.sites]
    SI <- sim$SI.count[target.gen, target.sites]
    A <- sim$A.count[target.gen, target.sites]
    AI <- sim$AI.count[target.gen, target.sites]
    N <- S+A
    prob <- matrix(c(S-SI, SI, A-AI, AI)/N, nrow=length(S))
    
    multisample <- {
        df <- as.data.frame(t(apply(prob, 1, function(p) rmultinom(1, size=nsample, prob=p))))
        names(df) <- c("SU", "SI", "AU", "AI")
        df
    }
    multisample <- multisample/nsample
    
    sampledf <- with(multisample,{
        transfun(data.frame(
            infected=SI+AI,
            sexual=SU+SI
        ))
    })
    
    sampledf
}

powerfun <- function(simlist,
                     level=0.05,
                     nsim=50,
                     nsample=50,
                     nsite=20,
                     transform=c("logit", "arcsin", "raw"),
                     test=test_spearman,
                     target.gen=1001, 
                     summarize=FALSE,
                     verbose=FALSE) {
    reslist <- vector("list", length(simlist))
    
    for (i in 1:length(simlist)) {
        if(verbose) cat(i, "\n")
        sim <- simlist[[i]]
        testlist <- vector("list", nsim)
        for (j in 1:nsim) {
            sample <- sample_sim(sim, nsample, nsite, transform, target.gen)
            testres <- suppressWarnings(test(sample)) 
            testlist[[j]] <- testres
        }
        comb_test <- do.call("rbind", testlist)
        rownames(comb_test) <- NULL
        
        reslist[[i]] <- comb_test
    }
    df <- do.call("rbind", reslist)
    
    if (summarize) {
        power <- sum(df[,2]<level, na.rm=TRUE)/(nsim*length(simlist))
        effect <- data.frame(
            median=median(df[,1], na.rm=TRUE),
            lwr=quantile(df[,1], 0.025, na.rm=TRUE),
            upr=quantile(df[,1], 0.975, na.rm=TRUE)
        )
        return(list(power=power,
             effect=effect))
    }
    return(df)
}
