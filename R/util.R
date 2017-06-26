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

outcross <- function(genotype, r) {
    recomb <- outer(c(-1, 1, 1, -1), c(genotype[1, 4], -genotype[2, 3]), "*")
    gamete <- colSums(genotype) + diag(genotype) + r * rowSums(recomb)
    new.genotype <- scaled_matrix(outer(gamete, gamete, "*"))
    if(scaled_sum(genotype) != 0) new.genotype <- new.genotype/scaled_sum(4*genotype)
    
    new.genotype
}

apply_gen <- function(x, name, FUN, 
                      gen=c(2001,2100)) {
    y <- x[name]
    lapply(y, function(z){
        if(is.vector(z)) {
            FUN(z[gen])
        } else if(is.matrix(z)) {
            FUN(z[gen,])
        } else if(is.array (z)) {
            FUN(z[gen,,])
        }
    })
}