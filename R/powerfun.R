find_power <- function(simlist,
                       s=0.5,
                       target.gen=1000,
                       samplesize=seq(20, 200, by=20),
                       level=0.05) {
    power <- matrix(0, nrow=dim(simlist)[2], ncol=length(samplesize))
    effsize <- matrix(0, nrow=dim(simlist)[2], ncol=length(samplesize))
    
    for(i in 1:dim(simlist)[2]) {
        multisample <- vector('list', dim(simlist)[1])
        testres <- vector('list', dim(simlist)[1])
        
        for (j in 1:dim(simlist)[1]) {
            sim <- simlist[j,i][[1]]
            SU <- sim$S.count[target.gen,]-sim$SI.count[target.gen,]
            SI <- sim$SI.count[target.gen,]
            AU <- sim$A.count[target.gen,]-sim$AI.count[target.gen,]
            AI <- sim$AI.count[target.gen,]
            
            prob <- matrix(c(s * SU, s*SI, s*SU, s*SI, AU, AI), nrow=length(SU))
            
            multisample[[j]] <- lapply(samplesize, function(s){
                apply(prob, 1, function(p) rmultinom(1, size=s, prob=p/sum(p)))
            })
            
            ## correlation between prevalence and male frequency
            testres[[j]] <- sapply(multisample[[j]], function(x) {
                cc <- cor.test(x[2,] + x[4,] + x[6,], x[1,])
                c(cor=cc$estimate, p.value=cc$p.value)
            })
        }
        power[i,] <- rowMeans(sapply(testres, function(x) x[2,]<level & !is.na(x[2,])))
        effmat <- sapply(testres, function(x) x[1,])
        effmat[is.na(effmat)] <- 0
        effsize[i,] <- rowMeans(effmat)
    }
    return(list(power=power, effsize=effsize))
}
