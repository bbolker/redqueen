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
