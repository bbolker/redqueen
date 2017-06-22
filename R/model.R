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
}

##' Lively2010
##' @param p allele frequencies at the first loci
##' @param q allele frequencies at second loci
lively_model <- function(p = c(0.1, 0.3, 0.6),
                         q = c(0.3, 0.3, 0.4),
                         s=0.5, 
                         r=0.2,
                         beta=5,
                         aU=0.001, aI=0.001,
                         bU=20, bI=3,
                         migrate.host=0.1, migrate.pathogen=0.02,
                         seed=NULL,
                         tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)
    
    S <- SI <- A <- AI <- array(0, dim=c(tmax, 3, 3))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax)
    
    migratefun <- function(p) {
        matrix(p, 3, 3)
    }
    
    S[1,,] <- 1000 * outer(p, q, "*")
    
    S.count[1] <- sum(S[1,,])
    N.count[1] <- sum(S[1,,])
    
    P <- 0
    
    for(t in 1:(tmax-1)){
        if (t==tburnin) {
            A[t,1,1] <- 1
            A.count[t] <- 1
        }
        
        S[t,,] <- S[t,,] + migratefun(migrate.host)
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * S[t,,] * (WI * P + WU * (1-P))
        p.gam <- colSums(Sp)/sum(Sp)
        q.gam <- rowSums(Sp)/sum(Sp)
        
        S[t+1,,] <- sum(Sp)*((1-r)*Sp/sum(Sp)+ r*outer(p.gam, q.gam, "*"))
        S.count[t+1] <- sum(S[t+1,,])
        
        A[t+1,,] <- A[t,,] * (WI * P + WU * (1-P))
        A.count[t+1] <- sum(A[t+1,,])
        AI[t,,] <- A[t,,] * P
        AI.count[t] <- sum(AI[t,,])
        
        SI[t,,] <- S[t,,] * P + migratefun(migrate.host)
        SI.count[t] <- sum(SI[t,,])
        
        N.count[t+1] <- S.count[t+1] + A.count[t+1]
        P <- 1 - exp(-beta*(SI[t,,]+AI[t,,])/N.count[t+1])
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, A=A, AI=AI, A.count=A.count, AI.count=AI.count
    ))
}


##' Discrete model without spatial structure
##' 
##' @param start list of starting values
##' @param s male offspring ratio
##' @param r host recombination probability
##' @param epsilon parasite mutation probability
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
discrete_model <- function(start=discrete_initialize(),
                           s=0.5, 
                           r=0.2, epsilon=0.01,
                           beta=5,
                           aU=0.001, aI=0.001,
                           bU=20, bI=3,
                           migrate.host=0.1, migrate.pathogen=0.02,
                           seed=NULL,
                           tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)

    S <- SI <- A <- AI <- array(0, dim=c(tmax, 4, 4))
    I <- matrix(0, nrow=tmax, ncol=4)
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax)
    
    migratefun <- function(p) {
        matrix(p, 4, 4)
    }
    
    S[1,,] <- start
    
    S.count[1] <- scaled_sum(S[1,,])
    N.count[1] <- scaled_sum(S[1,,])
    P <- 0

    for(t in 1:(tmax-1)){
        if (t==tburnin) {
            A[t,,] <- introduce()
            A.count[t] <- 1
        }
        
        S[t,,] <- S[t,,] + migratefun(migrate.host)
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * S[t,,] * (WI * P + WU * (1-P))
        S[t+1,,] <- outcross(Sp, r)
        S.count[t+1] <- scaled_sum(S[t+1,,])
        
        SI[t,,] <- S[t,,] * P
        SI.count[t] <- scaled_sum(SI[t,,])
        
        A[t+1,,] <- A[t,,] * (WI * P + WU * (1-P))
        A.count[t+1] <- scaled_sum(A[t+1,,])
        
        AI[t,,] <- A[t,,] * P
        AI.count[t] <- scaled_sum(AI[t,,])
        
        N.count[t+1] <- S.count[t+1] + A.count[t+1]
        
        I[t,] <- 1/2*(colSums(SI[t,,] +AI[t,,]) + diag(SI[t,,] +AI[t,,]))
        
        I[t,] <- ((1-epsilon)*I[t,]+epsilon * (sum(I[t,]) - (I[t,]+rev(I[t,])))/2) + migrate.pathogen
        
        inf <- scaled_matrix(outer(I[t,], I[t,], "+"))/2
        
        P <- 1 - exp(-beta*inf/N.count[t])
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, A=A, AI=AI, A.count=A.count, AI.count=AI.count, N.count=N.count, I=I
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

