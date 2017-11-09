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

stochastic_spatial_discrete_lim_model <- function(start,
                                                  n.site=40,
                                                  n.genotype=1,
                                                  epsilon.site=0.01,
                                                  c_b=1,
                                                  s=0.5, 
                                                  r.host=0.2,
                                                  beta=c(5, 5, 5, 5),
                                                  aU=0.001, aI=0.001,
                                                  bU=20, bI=3,
                                                  migrate.host=0.1, migrate.parasite=0.02,
                                                  r.parasite=0.05,
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
        SI[1,,,i] <- 0.04 * start
        S.count[1,i] <- scaled_sum(S[1,,,i])
        SI.count[1,i] <- scaled_sum(SI[1,,,i])
        N.count[1,i] <- scaled_sum(SI[1,,,i] + S[1,,,i])
    }
    
    inf <- array(0, dim=c(4, 4, n.site))
    P <- array(0, dim=c(4, 4, n.site))
    
    ratio <- array(0, dim=c(4, 4, n.site))
    
    asex.introduce <- function() {
        asex <- matrix(0, 4, 4)
        asex[which(.upr)[asex.genotype[sample(1:n.genotype, 1)]]] <- 1
        scaled_matrix(asex)
    }
    
    sex.migrate <- 1-(1-migrate.host)^(1/10)
    asex.migrate <- 1-(1-migrate.host * asex.mat)^(1/n.genotype)
    parasite.migrate <- 1-(1-migrate.parasite)^(1/4)
    
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
            
            Sp <- c_b * (1-s) *  ((S[t,,,i] - SI[t,,,i]) * WU +SI[t,,,i] * WI)
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
            I[t,,i] <- (1-r.parasite) * I.nomut + r.parasite/2 * (sum(I.nomut) - (I.nomut + rev(I.nomut))) + 
                as.numeric(runif(4) < parasite.migrate)
            
            lambda[t,,i] <- beta[i] * I[t,,i]
            
        }
        
        lambda.tot <- array(0, dim=c(4, n.site))
        ratio <- array(0, dim=c(4, 4, n.site))
        
        for(i in 1:n.site) {
            for(j in 1:n.site) {
                lambda.tot[,i] <- lambda.tot[,i] + ifelse(i==j, 1-epsilon.site, epsilon.site/(n.site-1)) * lambda[t,,j]
            }
            
            lambda.mat <- outer(lambda.tot[,i], lambda.tot[,i], "+")
            
            FOI <- lambda.mat/(2 * N.count[t+1,i])
            if (N.count[t+1,i] != 0) {
                P[,,i] <- 1 - exp(-FOI)
                ratio[,,i] <- lambda.tot[,i]/lambda.mat
                ratio[,,i][which(lambda.mat==0)] <- 0
                diag(ratio[,,i]) <- 1
                
                SI[t+1,,,i] <- binom_matrix(S[t+1,,,i], P[,,i])
                SI.count[t+1,i] <- scaled_sum(SI[t+1,,,i])
                
                AI[t+1,,,i] <- binom_matrix(A[t+1,,,i], P[,,i])
                AI.count[t+1,i] <- scaled_sum(AI[t+1,,,i])
            }
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
