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
                           aU=0.001, aI=0.001,
                           bU=20, bI=3,
                           migrate.host=0.1, migrate.pathogen=0.02,
                           epsilon=0.01,
                           seed=NULL,
                           tmax=1100, tburnin=500) {
    if(!is.null(seed)) set.seed(seed)
    
    S <- SI <- A <- AI <- P <- array(0, dim=c(tmax+1, 4, 4))
    I <- lambda <- matrix(0, nrow=tmax+1, ncol=4)
    N.count <- S.count <- SI.count <- A.count <- AI.count <- rep(0,tmax+1)
    
    S[1,,] <- start
    
    S.count[1] <- scaled_sum(S[1,,])
    N.count[1] <- scaled_sum(S[1,,])
    
    ratio <- 1
    
    for(t in 1:(tmax)){
        if (t==tburnin) {
            A[t,,] <- 10 * introduce()
            A.count[t] <- 10
        }
        
        WU <- bU/(1+aU*N.count[t])
        WI <- bI/(1+aI*N.count[t])
        
        Sp <- (1-s) * ((S[t,,]-SI[t,,])*WU + SI[t,,] * WI)
        S[t+1,,] <- pois_matrix(outcross(Sp, r)) + migratefun(migrate.host)
        S.count[t+1] <- scaled_sum(S[t+1,,])

        A[t+1,,] <- pois_matrix((A[t,,]- AI[t,,])*WU + AI[t,,] * WI)
        A.count[t+1] <- scaled_sum(A[t+1,,])
        
        N.count[t+1] <- S.count[t+1] + A.count[t+1]
        
        ## haplotype density
        tmp <- SI[t,,] + AI[t,,]
        
        I.nomut <- rowSums(tmp * ratio)
        I[t,] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + as.numeric(runif(4) <  migrate.pathogen)
        
        if(N.count[t+1] > 0) {
            lambda[t,] <- beta * I[t,]/(2*N.count[t+1])
        } else {
            lambda[t,] <- 0
        }
        
        FOI <- outer(lambda[t,], lambda[t,], "+")
        
        P[t+1,,] <- 1 - exp(-FOI)
        
        ratio <- lambda[t,]/FOI
        diag(ratio) <- 1
        ratio[which(is.nan(ratio))] <- 0
        
        SI[t+1,,] <- binom_matrix(S[t+1,,], P[t+1,,])
        SI.count[t+1] <- scaled_sum(SI[t+1,,])
        
        AI[t+1,,] <- binom_matrix(A[t+1,,], P[t+1,,])
        AI.count[t+1] <- scaled_sum(AI[t+1,,])
    }
    
    return(list(
        S=S, SI=SI, S.count=S.count, SI.count=SI.count, 
        A=A, AI=AI, A.count=A.count, AI.count=AI.count, 
        N.count=N.count, I=I, lambda=lambda,
        P=P
    ))
}

stochastic_spatial_discrete_model <- function(start,
                                   n.site=4,
                                   epsilon.site=0.01,
                                   s=0.5, 
                                   r.host=0.2,
                                   beta=c(5, 5, 5, 5),
                                   aU=0.001, aI=0.001,
                                   bU=20, bI=3,
                                   migrate.host=0.1, migrate.pathogen=0.02,
                                   epsilon=0.01,
                                   seed=NULL, simplify=TRUE,
                                   tmax=1100, tburnin=500) {
    if(!is.null(seed)) set.seed(seed)
    
    if(missing(start)) start <- discrete_initialize()
    
    S <- SI <- A <- AI <- array(0, dim=c(tmax+1, 4, 4, n.site))
    I <- lambda <- array(0, dim=c(tmax+1, 4, n.site))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- matrix(0, nrow=tmax+1, ncol=n.site)
    
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
                A[t,,,i] <- 10 * introduce()
                A.count[t,i] <- 10
            }
            
            WU <- bU/(1+aU*N.count[t,i])
            WI <- bI/(1+aI*N.count[t,i])
            
            Sp <- (1-s) *  ((S[t,,,i] - SI[t,,,i]) * WU +SI[t,,,i] * WI)
            S[t+1,,,i] <- pois_matrix(outcross(Sp, r.host)) + migratefun(migrate.host)
            S.count[t+1,i] <- scaled_sum(S[t+1,,,i])
            
            A[t+1,,,i] <- pois_matrix((A[t,,,i] - AI[t,,,i]) * WU + AI[t,,,i] * WI)
            A.count[t+1,i] <- scaled_sum(A[t+1,,,i])
            
            N.count[t+1,i] <- S.count[t+1,i] + A.count[t+1,i]
            
            tmp[,,i] <- SI[t,,,i] + AI[t,,,i]
            
            I.nomut <- rowSums(tmp[,,i] * ratio[,,i])
            I[t,,i] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + as.numeric(runif(4) <  migrate.pathogen)
            
            lambda[t,,i] <- beta[i] * I[t,,i]/(2 * N.count[t+1,i])
            lambda[t,,i][which(is.nan(lambda[t,,i]))] <- 0
            
        }
        
        lambda.tot <- array(0, dim=c(4, n.site))
        ratio <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            for(j in 1:n.site) {
                lambda.tot[,i] <- lambda.tot[,i] + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * lambda[t,,j]
            }
            
            FOI <- outer(lambda.tot[,i], lambda.tot[,i], "+")
            
            P[,,i] <- 1 - exp(-FOI)
            ratio[,,i] <- lambda.tot[,i]/FOI
            ratio[,,i][which(is.nan(ratio[,,i]))] <- 0
            diag(ratio[,,i]) <- 1

            SI[t+1,,,i] <- binom_matrix(S[t+1,,,i], P[,,i])
            SI.count[t+1,i] <- scaled_sum(SI[t+1,,,i])
            
            
            AI[t+1,,,i] <- binom_matrix(A[t+1,,,i], P[,,i])
            AI.count[t+1,i] <- scaled_sum(AI[t+1,,,i])
        }
        
    }
    
    if (simplify) {
        list(
            S.count=S.count, SI.count=SI.count,
            A.count=A.count, AI.count=AI.count
        )
    } else{
        list(
            S=S, SI=SI,
            S.count=S.count, SI.count=SI.count, 
            A=A, AI=AI,
            A.count=A.count, AI.count=AI.count, 
            N.count=N.count, I=I
        )
    }
}

stochastic_spatial_discrete_lim_model <- function(start,
                                                  n.site=4,
                                                  n.genotype=1,
                                                  epsilon.site=0.01,
                                                  s=0.5, 
                                                  r.host=0.1,
                                                  beta=c(5, 5, 5, 5),
                                                  aU=0.001, aI=0.001,
                                                  bU=20, bI=3,
                                                  migrate.host=0.1, migrate.pathogen=0.05,
                                                  epsilon=0.05,
                                                  seed=NULL, simplify=TRUE,
                                                  tmax=1100, tburnin=500) {
    if(!is.null(seed)) set.seed(seed)
    
    if(missing(start)) start <- discrete_initialize(N0=2000)
    
    S <- SI <- A <- AI <- array(0, dim=c(tmax+1, 4, 4, n.site))
    I <- lambda <- array(0, dim=c(tmax+1, 4, n.site))
    N.count <- S.count <- SI.count <- A.count <- AI.count <- matrix(0, nrow=tmax+1, ncol=n.site)
    
    asex.genotype <- sample(1:10, n.genotype)
    asex.mat <- matrix(0, 4, 4)
    asex.mat[which(.upr)[asex.genotype]] <- 1
    asex.mat <- scaled_matrix(asex.mat)
    
    for(i in 1:n.site) {
        S[1,,,i] <- start
        S.count[1,i] <- scaled_sum(S[1,,,i])
        N.count[1,i] <- scaled_sum(S[1,,,i])
    }
    
    inf <- array(0, dim=c(4, 4, n.site))
    P <- array(0, dim=c(4, 4, n.site))
    
    ratio <- array(0, dim=c(4, 4, n.site))
    
    asex.introduce <- function() {
        asex <- matrix(0, 4, 4)
        asex[which(.upr)[asex.genotype[sample(1:n.genotype, 1)]]] <- 1
        scaled_matrix(asex)
    }
    
    sex.migrate <- migrate.host
    asex.migrate <- migrate.host*asex.mat
    pathogen.migrate <- migrate.pathogen
    
    for(t in 1:tmax){
        
        tmp <- array(0, dim=c(4, 4, n.site))
        
        S.tmp <- array(0, dim=c(4, 4, n.site))
        A.tmp <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            if (t==tburnin) {
                A[t,,,i] <- 10 * asex.introduce()
                A.count[t,i] <- 10
                N.count[t,i] <- N.count[t,i] + 10
            }
            
            WU <- bU/(1+aU*N.count[t,i])
            WI <- bI/(1+aI*N.count[t,i])
            
            Sp <- (1-s) *  ((S[t,,,i] - SI[t,,,i]) * WU +SI[t,,,i] * WI)
            S.tmp[,,i] <- Sp
            
            A.tmp[,,i] <- (A[t,,,i] - AI[t,,,i]) * WU + AI[t,,,i] * WI
        }
        
        for(i in 1:n.site) {
            A.tmp2 <- S.tmp2 <- matrix(0, 4, 4)
            for(j in 1:n.site) {
                S.tmp2 <- S.tmp2 + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * S.tmp[,,j]
                A.tmp2 <- A.tmp2 + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * A.tmp[,,j]
            }
            
            S[t+1,,,i] <- pois_matrix(outcross(S.tmp2, r.host)) + migratefun(sex.migrate)
            S.count[t+1,i] <- scaled_sum(S[t+1,,,i])
            
            A[t+1,,,i] <- pois_matrix(A.tmp2)
            if (t >= tburnin) {
                A[t+1,,,i] <- A[t+1,,,i] + migratefun(asex.migrate)
            }
            
            
            A.count[t+1,i] <- scaled_sum(A[t+1,,,i])
            
            N.count[t+1,i] <- S.count[t+1,i] + A.count[t+1,i]
            
            tmp[,,i] <- SI[t,,,i] + AI[t,,,i]
            
            I.nomut <- rowSums(tmp[,,i] * ratio[,,i])
            I[t,,i] <- (1-epsilon) * I.nomut + epsilon/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + 
                as.numeric(runif(4) < pathogen.migrate)
            
            lambda[t,,i] <- beta[i] * I[t,,i]/(2 * N.count[t+1,i])
            lambda[t,,i][which(is.nan(lambda[t,,i]))] <- 0
        }
        
        lambda.tot <- array(0, dim=c(4, n.site))
        ratio <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            for(j in 1:n.site) {
                lambda.tot[,i] <- lambda.tot[,i] + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * lambda[t,,j]
            }
            
            FOI <- outer(lambda.tot[,i], lambda.tot[,i], "+")
            
            P[,,i] <- 1 - exp(-FOI)
            ratio[,,i] <- lambda.tot[,i]/FOI
            ratio[,,i][which(is.nan(ratio[,,i]))] <- 0
            diag(ratio[,,i]) <- 1
            
            SI[t+1,,,i] <- binom_matrix(S[t+1,,,i], P[,,i])
            SI.count[t+1,i] <- scaled_sum(SI[t+1,,,i])
            
            AI[t+1,,,i] <- binom_matrix(A[t+1,,,i], P[,,i])
            AI.count[t+1,i] <- scaled_sum(AI[t+1,,,i])
        }
        
    }
    
    if (simplify) {
        list(
            S.count=S.count, SI.count=SI.count,
            A.count=A.count, AI.count=AI.count
        )
    } else{
        list(
            S=S, SI=SI,
            S.count=S.count, SI.count=SI.count, 
            A=A, AI=AI,
            A.count=A.count, AI.count=AI.count, 
            N.count=N.count, I=I
        )
    }
}
