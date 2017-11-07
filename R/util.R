##' Scale matrix
scaled_matrix <- function(x) {
    x <- x + t(x)
    diag(x) <- diag(x)/2
    x
}

.upr <- upper.tri(matrix(0, 4, 4), diag=TRUE)

##' Scale 
scaled_sum <- function(x) {
    sum(x[.upr])
}

introduce <- function() {
    new.genotype <- scaled_matrix(matrix(rmultinom(1, size=1, prob=.upr), 4, 4))
    new.genotype
}

##' Initialize discrete model
##' 
##' @param p proportion of 
##' @param N0 initial population size
##' @param I0 initial number of infected
discrete_initialize <- function(p=0.5, q=0.5, 
                                N0=2000) {
    locus1 <- c(p, 1-p)
    locus2 <- c(q, 1-q)
    
    gamete <- unlist(lapply(locus1, function(x) x*locus2))
    genotype <- outer(gamete, gamete, "*")
    genotype <- scaled_matrix(genotype)
    
    S <- N0 * genotype
    return(S)
}

outcross <- function(genotype, r) {
    recomb <- outer(c(-1, 1, 1, -1), c(genotype[1, 4], -genotype[2, 3]), "*")
    gamete <- colSums(genotype) + diag(genotype) + r * rowSums(recomb)
    new.genotype <- scaled_matrix(outer(gamete, gamete, "*"))
    if(scaled_sum(genotype) != 0) new.genotype <- new.genotype/scaled_sum(4*genotype)
    
    new.genotype
}

