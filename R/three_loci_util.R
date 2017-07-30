##' Scale matrix
scaled_matrix <- function(x) {
    x <- x + t(x)
    diag(x) <- diag(x)/2
    x
}

##' Scale 
scaled_sum <- function(x) {
    sum(x[.upr])
}

.upr <- upper.tri(matrix(0, 8, 8), diag=TRUE)

introduce3 <- function() {
    new.genotype <- scaled_matrix(matrix(rmultinom(1, size=1, prob=.upr), 8, 8))
    new.genotype
}

##' Initialize discrete model
##' 
##' @param p proportion of 
##' @param N0 initial population size
##' @param I0 initial number of infected
discrete_initialize3 <- function(p1=0.5, p2=0.5, 
                                p3=0.5,
                                N0=2000) {
    locus1 <- c(p1, 1-p1)
    locus2 <- c(p2, 1-p2)
    locus3 <- c(p3, 1-p3)
    ## see Hodgson & Otto (2012)
    gamete <- unlist(lapply(locus3, "*", unlist(lapply(locus1, "*", locus2))))
    genotype <- outer(gamete, gamete, "*")
    genotype <- scaled_matrix(genotype)
    
    S <- N0 * genotype
    return(S)
}

## taken from Hodgson and Otto (2012)
recombination <- function(xh, r) {
    ##Indices
    ir.2a <- array(c(-1,1,1,-1,0,0,0,0))
    ir.2b <- array(c(0,0,0,0,-1,1,1,-1))
    
    ir.3a <- array(c(0,1,0,-1,0,-1,0,1))
    ir.3b <- array(c(-1,0,1,0,1,0,-1,0))
    
    ir.4a <- array(c(-1,1,0,0,1,-1,0,0))
    ir.4b <- array(c(0,0,1,-1,0,0,-1,1))
    
    ir.5a <- array(c(0,1,-1,0,0,-1,1,0))
    ir.5b <- array(c(-1,0,0,1,1,0,0,-1))
    
    ir.6a <- array(c(1,0,-1,0,0,-1,0,1))
    ir.6b <- array(c(0,1,0,-1,-1,0,1,0))
    
    ir.7a <- array(c(1,0,0,-1,-1,0,0,1))
    ir.7b <- array(c(0,-1,1,0,0,1,-1,0))
    
    D <- function(x1, x2, x3, x4) {
        xh[x1]*xh[x2] - xh[x3]*xh[x4]
    }	
    
    xh <- (xh 	+ ir.2a*r*D(1,4,2,3) ## A B (C)
           + ir.2b*r*D(5,8,6,7) ## A B (c)
           + ir.3a*(2*r - 2*r^2)*D(4,6,2,8) ## A B + B C - both
           + ir.3b*(2*r - 2*r^2)*D(1,7,3,5) ## A B + B C - both
           + ir.4a*r*D(1,6,2,5) ## B C (same A)
           + ir.4b*r*D(4,7,3,8) ## B C
           + ir.5a*(r - r^2)*D(3,6,2,7) ## B C - both
           + ir.5b*(r - r^2)*D(1,8,4,5) ## B C - both
           + ir.6a*r*D(3,6,1,8) ## A B
           + ir.6b*r*D(4,5,2,7) ## A B
           + ir.7a*r^2*D(2,7,3,6) ## both
           + ir.7b*r^2*D(4,5,1,8)) ## both
    xh
}

outcross3 <- function(genotype, r) {
    gamete.norecomb <- (colSums(genotype)+diag(genotype))/2
    if(scaled_sum(genotype) != 0) gamete.norecomb <- gamete.norecomb/scaled_sum(genotype)
    
    gamete <- recombination(gamete.norecomb, r)
    new.genotype <- scaled_matrix(outer(gamete, gamete, "*")) * scaled_sum(genotype)

    new.genotype
}

pathogen.mutate <- function(haplotype, epsilon) {
    ## h1 <- haplotype[1:4]
    ## h2 <- haplotype[5:8]
    ## (1-epsilon) * haplotype + epsilon*c(sum(h1)-(h1+rev(h1))+h2, sum(h2)-(h2+rev(h2))+h1)/3
    
    ## this is faster
    (1-epsilon) * haplotype + epsilon * c(
        haplotype[2]+haplotype[3]+haplotype[5],
        haplotype[1]+haplotype[4]+haplotype[6],
        haplotype[1]+haplotype[4]+haplotype[7],
        haplotype[2]+haplotype[3]+haplotype[8],
        haplotype[6]+haplotype[7]+haplotype[1],
        haplotype[5]+haplotype[8]+haplotype[2],
        haplotype[5]+haplotype[8]+haplotype[3],
        haplotype[6]+haplotype[7]+haplotype[4]
    )/3
}
