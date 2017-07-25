lm2 <- function(sim, 
                nsample,
                nsite,
                target.gen) {
    sites <- ncol(sim$S.count)
    target.sites <- sample(1:sites, nsite)
    
    S <- sim$S.count[target.gen, target.sites]
    SI <- sim$SI.count[target.gen, target.sites]
    A <- sim$A.count[target.gen, target.sites]
    AI <- sim$AI.count[target.gen, target.sites]
    N <- S+A
    prob <- matrix(c(S-SI, SI, A-AI, AI)/N, nrow=length(S))
    
    multisample <- lapply(nsample, function(s){
        df <- as.data.frame(t(apply(prob, 1, function(p) rmultinom(1, size=s, prob=p))))
        names(df) <- c("SU", "SI", "AU", "AI")
        df
    })
    
    testres <- sapply(multisample, function(x) {
        cc <- cor.test((x$SI + x$AI)/nsample, (x$SU + x$SI)/nsample)
        c(cor=cc$estimate, p.value=cc$p.value)
    })
    
    as.data.frame(t(testres))
}

powerfun <- function(simlist,
                     level=0.05,
                     nsim=50,
                     nsample=50,
                     nsite=20,
                     target.gen=1001) {
    reslist <- vector("list", length(simlist))
    
    for(i in 1:length(simlist)) {
        sim <- simlist[[i]]
        reslist[[i]] <- do.call("rbind", replicate(nsim, lm2(sim, nsample, nsite, target.gen), simplify=FALSE))
    }
    df <- do.call("rbind", reslist)
    
    power <- mean(df[,2]<level)
    effect <- data.frame(
        median=median(df[,1]),
        lwr=quantile(df[,1], 0.025),
        upr=quantile(df[,1], 0.975)
    )
    list(power=power,
         effect=effect)
}
