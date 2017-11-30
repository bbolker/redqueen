test_lm <- function(sample) {
    fit <- cor.test(sample$infected, sample$sexual, method="pearson")
    data.frame(
        effect.size=fit$estimate,
        p.value=fit$p.value
    )
}

test_quad <- function(sample) {
    fit <- lm(sexual~I(infected)+I(infected^2), data=sample)
    ss <- summary(fit)
    data.frame(
        effect.size=fit$coefficients[3],
        p.value=ifelse(is.na(fit$coefficients[3]), NA, summary(fit)$coefficients[3,'Pr(>|t|)']) 
    )
}

test_quad_rq <- function(sample, tau=0.9) {
    fit <- try(quantreg::rq(sexual~I(infected)+I(infected^2), data=sample, tau=tau), silent=TRUE)
    if (inherits(fit, "try-error") || fit$coefficients[3]==0) {
        data.frame(
            effect.size=NA,
            p.value=NA
        )
    } else {
        ss <- try(summary(fit, se="boot"), silent=TRUE)
        
        data.frame(
            effect.size=ifelse(inherits(ss, "try-error"), NA, ss$coefficients[3,3]),
            p.value=ifelse(inherits(ss, "try-error"), NA, ss$coefficients[3,4])
        )
    }
}

test_spearman <- function(sample) {
    fit <- cor.test(sample$infected, sample$sexual, method="spearman")
    data.frame(
        effect.size=fit$estimate,
        p.value=fit$p.value
    )
}

.logitfun <- function(p) {
    up <- unlist(p)
    minp <- min(c(up[which(up != 0)], 1-max(up[which(up != 1)])))
    
    if(any(p==0) || any(p==1)) {
        log((p+minp)/(1-p+minp))
    } else {
        log(p/(1-p))
    }
} 

transfun <- function(transform=c("logit", "arcsin", "raw")) {
    transform <- match.arg(transform)
    switch(transform,
        logit={
            function(df) {
                as.data.frame(lapply(df, .logitfun))
            } 
        },
        arcsin={
            function(df) asin(sqrt(df))
        },
        raw={
            function(df) df
        }
    )
}

sample_sim <- function(sim,
                  nsample,
                  nsite,
                  transform=c("raw", "arcsin", "logit"),
                  target.gen,
                  target.sites) {
    transfun <- transfun(transform)
    
    sites <- ncol(sim$S.count)
    
    if(missing(target.sites)) target.sites <- sample(1:sites, nsite)
    
    S <- sim$S.count[target.gen, target.sites]
    SI <- sim$SI.count[target.gen, target.sites]
    A <- sim$A.count[target.gen, target.sites]
    AI <- sim$AI.count[target.gen, target.sites]
    N <- S+A
    
    if (length(target.gen) == 1) {
        prob <- matrix(
            c(S-SI,SI,A-AI,AI)/N, 
            nrow=length(S)
        )
    } else {
        prob <- matrix(
            c(
                colMeans((S-SI)/N),
                colMeans(SI/N),
                colMeans((A-AI)/N),
                colMeans(AI/N)
            ), 
            nrow=ncol(S)
        )
    }

    multisample <- {
        df <- as.data.frame(t(apply(prob, 1, function(p) rmultinom(1, size=nsample, prob=p))))
        names(df) <- c("SU", "SI", "AU", "AI")
        df
    }
    multisample <- multisample/nsample
    
    sampledf <- transfun(with(multisample,{
        data.frame(
            infected=SI+AI,
            sexual=SU+SI
        )
    }))
    
    sampledf
}

powerfun <- function(simlist,
                     level=0.05,
                     nsim=50,
                     nsample=50,
                     nsite=20,
                     transform=c("raw", "arcsin", "logit"),
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

