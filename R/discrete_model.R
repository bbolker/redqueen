##' Initialize discrete model
##' 
##' @param p proportion of 
##' @param N0 initial population size
##' @param I0 initial number of infected
discrete_initialize <- function(p=0.5, q=0.5, 
                                N0=1000, I0=1) {
    locus1 <- c(p, 1-p)
    locus2 <- c(q, 1-q)
    
    gamete <- unlist(lapply(locus1, function(x) x*locus2))
    genotype <- outer(gamete, gamete, "*")
    genotype <- scaled_matrix(genotype)
    
    S <- (N0-I0) * genotype
    SI <- I0 * genotype
    
    list(S=S, SI=SI)
}

##' Discrete model without spatial structure
##' 
##' @param start list of starting values
##' @param s male offspring ratio
##' @param r.host host recombination probability
##' @param r.pathogen pathogen recombination probability
##' @param beta realized fecundity of parasite (Lively, 2010)
##' @param aU scales the effect of total host density on offpspring production by uninfected hosts
##' @param aI scales the effect of total host density on offpspring production by infected hosts
##' @param bU number of offspring produced by uninfected hosts in the absence of competitors
##' @param bI number of offspring produced by infected hosts in the absence of competitors
##' @param migrate.host sexual host migration probability per genotype
##' @param migrate.pathogen pathogen migration probability per genotype
##' @param tmax maximum number of generations
##' @param tburnin burn-in period
##' @return a list containing simulation results
discrete_model <- function(start,
                           s=0.5, 
                           r.host=0.2, r.pathogen=0.1,
                           beta=20,
                           aU=0.001, aI=0.001,
                           bU=20, bI=1,
                           migrate=c("none","deterministic","stochastic"),
                           migrate.host=0.1, migrate.pathogen=0.02,
                           seed=NULL,
                           tmax=2100, tburnin=1000) {
    migrate <- match.arg(migrate)
    
    if(!is.null(seed)) set.seed(seed)

    if(missing(start)) start <- discrete_initialize()
    
    S <- SI <- A <- AI <- I <- array(0, dim=c(tmax, 4, 4))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- I.count <- rep(0,tmax)
    
    outcross <- function(genotype, r) {
        recomb <- outer(c(-1, 1, 1, -1), c(genotype[1, 4], -genotype[2, 3]), "*")
        gamete <- colSums(genotype) + diag(genotype) + r * rowSums(recomb)
        new.genotype <- scaled_matrix(outer(gamete, gamete, "*"))
        if(scaled_sum(genotype) != 0) new.genotype <- new.genotype/scaled_sum(4*genotype)
        
        new.genotype
    }
    
    migratefun <- function(p) {
        switch(migrate,
               none=matrix(0, 4, 4),
               deterministic=matrix(p, 4, 4),
               stochastic=scaled_matrix(matrix(as.numeric(runif(16)<p*.upr), 4, 4)))
    }
    
    S[1,,] <- start$S
    SI[1,,] <- start$SI
    
    S.count[1] <- scaled_sum(S[1,,])
    SI.count[1] <- scaled_sum(SI[1,,])

    for(t in 1:(tmax-1)){
        if (t==tburnin) {
            A[t,,] <- introduce()
            A.count[t] <- 1
        }
        
        
        S[t,,] <- S[t,,] + migratefun(migrate.host)
        SI[t,,] <- SI[t,,] + migratefun(migrate.pathogen)
        
        I[t,,] <- outcross(SI[t,,] +AI[t,,], r.pathogen)
        I.count[t] <- scaled_sum(I[t,,])
        
        N.count[t] <- S.count[t] + A.count[t] + I.count[t]
        
        ## rho <- (SI[t,,] + AI[t,,])
        ## if(I.count[t] != 0){
        ##    rho <- rho/I.count[t]
        ## }
        ## rho.recomb <- outer(c(-1, 1, 1, -1), c(rho[1, 4], -rho[2,3]), "*")
        ## rho.gamete <- colSums(rho) + diag(rho) + r * rowSums(rho.recomb)
        ## rho.gamete <- rho.gamete/2
        ## scaled_matrix(outer(rho.gamete, rho.gamete, "*"))    
        
        P <- 1 - exp(-beta*I[t,,]/N.count[t])
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * S[t,,] * (WI * P + WU * (1-P))
        S.count[t+1] <- scaled_sum(Sp)
        S[t+1,,] <- outcross(Sp, r.host)
        
        ## recomb.gamete <- outer(c(-1, 1, 1, -1), c(Sp[1,4],-Sp[2,3]), "*")
        ## S.gamete <- colSums(Sp) + diag(Sp) + r * rowSums(recomb.gamete)
        ## if(sum(S.gamete)!=0) {
        ##     S.gamete <- S.gamete/sum(S.gamete)
        ## }
        ## S.outcross <- scaled_matrix(outer(S.gamete, S.gamete, "*"))
        ## S[t+1,,] <- S.count[t+1] * S.outcross
        
        A[t+1,,] <- A[t,,] * (WI * P + WU * (1-P))
        A.count[t+1] <- scaled_sum(A[t+1,,])
        AI[t+1,,] <- A[t,,] * P
        AI.count[t+1] <- scaled_sum(AI[t+1,,])
        
        SI[t+1,,] <- S[t,,] * P
        SI.count[t+1] <- scaled_sum(SI[t+1,,])
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, A=A, AI=AI, A.count=A.count, AI.count=AI.count, I=I, I.count=I.count
    ))
}

spatial_discrete_model <- function(start,
                                   n.site=10, epsilon=0.001,
                                   s=0.5, 
                                   r.host=0.2, r.pathogen=0.1,
                                   beta=20,
                                   aU=0.001, aI=0.001,
                                   bU=20, bI=1,
                                   migrate=c("none","deterministic","stochastic"),
                                   migrate.host=0.1, migrate.pathogen=0.02,
                                   seed=NULL,
                                   tmax=2100, tburnin=1000) {
    migrate <- match.arg(migrate)
    
    if(!is.null(seed)) set.seed(seed)
    
    if(missing(start)) start <- discrete_initialize()
    
    S <- SI <- A <- AI <- I <- array(0, dim=c(tmax, 4, 4, n.site))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- I.count <- matrix(0, nrow=tmax, ncol=n.site)
    
    introduce <- function() {
        new.genotype <- scaled_matrix(matrix(rmultinom(1, size=1, prob=.upr), 4, 4))
        new.genotype
    }
    
    migratefun <- function(p) {
        switch(migrate,
               none=matrix(0, 4, 4),
               deterministic=matrix(p, 4, 4),
               stochastic=scaled_matrix(matrix(as.numeric(runif(16)<p*.upr), 4, 4)))
    }
    
    for(i in 1:n.site) {
        S[1,,,i] <- start$S
        SI[1,,,i] <- start$SI
        S.count[1,i] <- scaled_sum(S[1,,,i])
        SI.count[1,i] <- scaled_sum(SI[1,,,i])
    }
    
    for(t in 1:(tmax-1)){
        ## calculate infection stuff first...
        for(i in 1:n.site) {
            if (t==tburnin) {
                A[t,,,i] <- introduce()
                A.count[t,i] <- 1
            }
            
            
            S[t,,,i] <- S[t,,,i] + migratefun(migrate.host)
            SI[t,,,i] <- SI[t,,,i] + migratefun(migrate.pathogen)
            
            I[t,,,i] <- outcross(SI[t,,,i] +AI[t,,,i], r.pathogen)
            I.count[t,i] <- scaled_sum(I[t,,,i])
            
            N.count[t,i] <- S.count[t,i] + A.count[t,i] + I.count[t,i]
        }
        
        
        
        for(i in 1:n.site) {
            I.tot <- matrix(0, 4, 4)
            
            for(j in 1:n.site) {
                I.tot <- I.tot + ifelse(i==j, 1, epsilon) * I[t,,,j]
            }
            
            P <- 1 - exp(-beta*I.tot/N.count[t,i])
            
            WU <- bU/(1+aU*N.count[t,i])
            WI <- bI/(1+aI*N.count[t,i])
            
            Sp <- (1-s) * S[t,,,i] * (WI * P + WU * (1-P))
            S.count[t+1,i] <- scaled_sum(Sp)
            S[t+1,,,i] <- outcross(Sp, r.host)
            
            A[t+1,,,i] <- A[t,,,i] * (WI * P + WU * (1-P))
            A.count[t+1,i] <- scaled_sum(A[t+1,,,i])
            AI[t+1,,,i] <- A[t,,,i] * P
            AI.count[t+1,i] <- scaled_sum(AI[t+1,,,i])
            
            SI[t+1,,,i] <- S[t,,,i] * P
            SI.count[t+1,i] <- scaled_sum(SI[t+1,,,i])
        }
        
        
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, A=A, AI=AI, A.count=A.count, AI.count=AI.count, I=I, I.count=I.count
    ))
}

