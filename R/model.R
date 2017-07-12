##' Lively2010
##' @param p allele frequencies at the first loci
##' @param q allele frequencies at second loci
lively_model <- function(p = c(0.33, 0.33, 0.34),
                         q = c(0.33, 0.33, 0.34),
                         s=0.5, 
                         r=0.2,
                         beta=5,
                         aU=0.001, aI=0.001,
                         bU=20, bI=3,
                         migrate.host=0.1, migrate.pathogen=0.02,
                         seed=NULL,
                         tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)
    
    S <- SI <- A <- AI <- array(0, dim=c(tmax+1, 3, 3))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax+1)
    
    migratefun <- function(p) {
        matrix(as.numeric(runif(9)<p), 3, 3)
    }
    
    S[1,,] <- 1000 * outer(p, q, "*")
    
    S.count[1] <- sum(S[1,,])
    N.count[1] <- sum(S[1,,])
    
    P <- 1e-4
    
    for(t in 1:(tmax)){
        if (t==tburnin) {
            A[t,1,1] <- 1
            A.count[t] <- 1
        }
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * S[t,,] * (WI * P + WU * (1-P))
        p.gam <- colSums(Sp)/sum(Sp)
        q.gam <- rowSums(Sp)/sum(Sp)
        
        S[t+1,,] <- sum(Sp)*((1-r)*Sp/sum(Sp)+ r*outer(p.gam, q.gam, "*")) + migratefun(migrate.host)
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
                           r=0.2,
                           beta=5,
                           aU=0.001, aI=0.001,
                           bU=20, bI=3,
                           migrate.host=0.1, migrate.pathogen=0.02,
                           epsilon=0.01,
                           seed=NULL,
                           tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)

    S <- SI <- A <- AI <- array(0, dim=c(tmax+1, 4, 4))
    I <- lambda <- matrix(0, nrow=tmax+1, ncol=4)
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax+1)
    
    migratefun <- function(p) {
        matrix(p, 4, 4)
    }
    
    S[1,,] <- start
    
    S.count[1] <- scaled_sum(S[1,,])
    N.count[1] <- scaled_sum(S[1,,])
    P <- 0
    ratio <- 0

    for(t in 1:(tmax)){
        if (t==tburnin) {
            A[t,,] <- introduce()
            A.count[t] <- 1
        }
        
        S[t,,] <- S[t,,]
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * S[t,,] * (WI * P + WU * (1-P)) 
        S[t+1,,] <- outcross(Sp, r) + migratefun(migrate.host)
        S.count[t+1] <- scaled_sum(S[t+1,,])
        
        SI[t,,] <- S[t,,] * P
        SI.count[t] <- scaled_sum(SI[t,,])
        
        A[t+1,,] <- A[t,,] * (WI * P + WU * (1-P))
        A.count[t+1] <- scaled_sum(A[t+1,,])
        
        AI[t,,] <- A[t,,] * P
        AI.count[t] <- scaled_sum(AI[t,,])
        
        N.count[t+1] <- S.count[t+1] + A.count[t+1]
        
        ## haplotype density
        tmp <- SI[t,,] + AI[t,,]
        
        I.nomut <- rowSums(tmp * ratio)
        I[t,] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + migrate.pathogen
        
        lambda[t,] <- beta * I[t,]/(2*N.count[t+1])
        
        inf <- outer(lambda[t,], lambda[t,], "+")
        
        P <- 1 - exp(-inf)
        
        ratio <- lambda[t,]/inf
        diag(ratio) <- 1
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, 
        A=A, AI=AI, A.count=A.count, AI.count=AI.count, 
        N.count=N.count, I=I, lambda=lambda
    ))
}

spatial_discrete_model <- function(start,
                                   n.site=10,
                                   epsilon.site=0.01, epsilon=0.01,
                                   s=0.5,
                                   r.host=0.2,
                                   beta=5,
                                   aU=0.001, aI=0.001,
                                   bU=20, bI=3,
                                   migrate.host=0.1, migrate.pathogen=0.02,
                                   seed=NULL,
                                   tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)
    
    if(missing(start)) start <- discrete_initialize()
    
    S <- SI <- A <- AI <- array(0, dim=c(tmax+1, 4, 4, n.site))
    I <- lambda <- array(0, dim=c(tmax+1, 4, n.site))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- matrix(0, nrow=tmax+1, ncol=n.site)
    
    migratefun <- function(p) {
        matrix(p, 4, 4)
    }
    
    for(i in 1:n.site) {
        S[1,,,i] <- start
        S.count[1,i] <- scaled_sum(S[1,,,i])
        N.count[1,i] <- scaled_sum(S[1,,,i])
    }
    
    inf <- array(0, dim=c(4, 4, n.site))
    P <- array(0, dim=c(4, 4, n.site))
    
    ratio <- array(0, dim=c(4, 4, n.site))
    
    for(t in 1:(tmax)){
        
        tmp <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            if (t==tburnin) {
                A[t,,,i] <- introduce()
                A.count[t,i] <- 1
            }
            
            WU <- bU/(1+aU*N.count[t,i])
            WI <- bI/(1+aI*N.count[t,i])
            
            Sp <- (1-s) * S[t,,,i] * (WI * P[,,i] + WU * (1-P[,,i]))
            S[t+1,,,i] <- outcross(Sp, r.host) + migratefun(migrate.host)
            S.count[t+1,i] <- scaled_sum(S[t+1,,,i])
            
            SI[t,,,i] <- S[t,,,i] * P[,,i]
            SI.count[t,i] <- scaled_sum(SI[t,,,i])
            
            A[t+1,,,i] <- A[t,,,i] * (WI * P[,,i] + WU * (1-P[,,i]))
            A.count[t+1,i] <- scaled_sum(A[t+1,,,i])
            
            AI[t,,,i] <- A[t,,,i] * P[,,i]
            AI.count[t,i] <- scaled_sum(AI[t,,,i])
            
            N.count[t+1,i] <- S.count[t+1,i] + A.count[t+1,i]
            
            tmp[,,i] <- SI[t,,,i] + AI[t,,,i]
            
            I.nomut <- rowSums(tmp[,,i] * ratio[,,i])
            
            I[t,,i] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + migrate.pathogen
            
            lambda[t,,i] <- beta[i] * I[t,,i]/(2*N.count[t+1,i])
            
        }
        
        lambda.tot <- array(0, dim=c(4, n.site))
        ratio <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            for(j in 1:n.site) {
                lambda.tot[,i] <- lambda.tot[,i] + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * lambda[t,,j]
            }
            
            inf <- outer(lambda.tot[,i], lambda.tot[,i], "+")
            
            P[,,i] <- 1 - exp(-inf)
            ratio[,,i] <- lambda[t,,i]/inf
            diag(ratio[,,i]) <- 1
        }
        
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, A=A, AI=AI, A.count=A.count, AI.count=AI.count, N.count=N.count, I=I
    ))
}
