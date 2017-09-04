test_lm <- function(sample) {
    fit <- lm(sexual~infected, data=sample)
    data.frame(
        effect.size=fit$coefficients[2],
        p.value=ifelse(is.na(fit$coefficients[2]), NA, summary(fit)$coefficients[2,'Pr(>|t|)']) 
    )
}

test_quad <- function(sample) {
    fit <- lm(sexual~I(infected)+I(infected^2), data=sample)
    data.frame(
        effect.size=fit$coefficients[3],
        p.value=ifelse(is.na(fit$coefficients[3]), NA, summary(fit)$coefficients[3,'Pr(>|t|)']) 
    )
}

test_quad_rq <- function(sample, tau=0.9) {
    fit <- quantreg::rq(sexual~I(infected)+I(infected^2), data=sample, tau=tau)
    data.frame(
        effect.size=fit$coefficients[3],
        p.value=NA
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
            function(p) asin(sqrt(p))
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
                     test=list(spearman=test_spearman),
                     target.gen=1001,
                     verbose=FALSE) {
    reslist <- vector("list", length(simlist))
    
    for (i in 1:length(simlist)) {
        if(verbose) cat(i, "\n")
        sim <- simlist[[i]]
        testlist <- vector("list", nsim)
        for (j in 1:nsim) {
            sample <- sample_sim(sim, nsample, nsite, transform, target.gen)
            suppressWarnings(testres <- lapply(test, function(test) test(sample)))
            tt <- do.call("rbind", testres)
            tt$test <- names(testres)
            testlist[[j]] <- tt
        }
        comb_test <- do.call("rbind", testlist)
        rownames(comb_test) <- NULL
         
        reslist[[i]] <- comb_test
    }
    df <- do.call("rbind", reslist)
    
    return(df)
}

