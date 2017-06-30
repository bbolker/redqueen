##' Discrete stochsatic model without spatial structure
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
##' @param nb.size size parameter for negative binomial
##' @param tmax maximum number of generations
##' @param tburnin burn-in period
##' @return a list containing simulation results
stochastic_discrete_model <- function(start=discrete_initialize(),
                           s=0.5, 
                           r=0.2,
                           beta=5,
                           aU=0.0001, aI=0.0001,
                           bU=30, bI=3,
                           migrate.host=0.1, migrate.pathogen=0.02,
                           epsilon=0.01,
                           seed=NULL,
                           tmax=2100, tburnin=1000) {
    if(!is.null(seed)) set.seed(seed)
    
    S <- SI <- A <- AI <- P <- array(0, dim=c(tmax+1, 4, 4))
    I <- lambda <- matrix(0, nrow=tmax+1, ncol=4)
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax+1)
    
    migratefun <- function(p) {
        scaled_matrix(matrix(runif(16) < c(p * .upr), 4, 4))
    }
    
    pois_matrix <- function(mat){
        tmpmat <- .upr
        rvec <- rpois(n=10, lambda=mat[.upr])
        tmpmat[.upr] <- rvec
        
        scaled_matrix(tmpmat)
    }
    
    binom_matrix <- function(mat, p) {
        tmpmat <- .upr
        rvec <- rbinom(n=10, size=mat[.upr], prob=p[.upr])
        tmpmat[.upr] <- rvec
        
        scaled_matrix(tmpmat)
    }
    
    S[1,,] <- start
    
    S.count[1] <- scaled_sum(S[1,,])
    N.count[1] <- scaled_sum(S[1,,])
    
    SI.zero <- scaled_matrix(matrix(1, 4, 4))
    
    I.zero <- c(4, 4, 4, 4)
    
    lambda.zero <- beta * I.zero/(2*N.count[1])
    
    FOI <- outer(lambda.zero, lambda.zero, "+")
    
    P[1,,] <- 1 - exp(-FOI)
    
    ratio <- lambda.zero/FOI
    diag(ratio) <- 1
    
    for(t in 1:(tmax)){
        if (t==tburnin) {
            A[t,,] <- introduce()
            A.count[t] <- 1
        }
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        SI[t,,] <- binom_matrix(S[t,,], P[t,,])
        SI.count[t] <- scaled_sum(SI[t,,])
        
        Sp <- (1-s) * ((S[t,,]-SI[t,,])*WU + SI[t,,] * WI)
        S[t+1,,] <- pois_matrix(outcross(Sp, r)) + migratefun(migrate.host)
        S.count[t+1] <- scaled_sum(S[t+1,,])
        
        AI[t,,] <- binom_matrix(A[t,,], P[t,,])
        AI.count[t] <- scaled_sum(AI[t,,])
        
        A[t+1,,] <- pois_matrix((A[t,,]- AI[t,,])*WU + AI[t,,] * WI)
        A.count[t+1] <- scaled_sum(A[t+1,,])
        
        N.count[t+1] <- S.count[t+1] + A.count[t+1]
        
        ## haplotype density
        tmp <- SI[t,,] + AI[t,,]
        
        I.nomut <- rowSums(tmp * ratio)
        I[t,] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + as.numeric(runif(4) < 0.1)
        
        lambda[t,] <- beta * I[t,]/(2*N.count[t+1])
        
        FOI <- outer(lambda[t,], lambda[t,], "+")
        
        P[t+1,,] <- 1 - exp(-FOI)
        
        ratio <- lambda[t,]/FOI
        diag(ratio) <- 1
        ratio[which(is.nan(ratio))] <- 0
        
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, 
        A=A, AI=AI, A.count=A.count, AI.count=AI.count, 
        N.count=N.count, I=I, lambda=lambda,
        P=P
    ))
}
