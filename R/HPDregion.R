## https://stackoverflow.com/questions/24198514/ggplot2-modify-geom-density2d-to-accept-weights-as-a-parameter
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
    nx <- length(x)
    if (length(y) != nx) 
        stop("data vectors must be the same length")
    if (length(w) != nx && length(w) != 1)
        stop("weight vectors must be 1 or length of data")
    gx <- seq(lims[1], lims[2], length = n) # gridpoints x
    gy <- seq(lims[3], lims[4], length = n) # gridpoints y
    if (missing(h)) 
        h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
    if (missing(w)) 
        w <- numeric(nx)+1;
    h <- h/4
    ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
    ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
    z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
    return(list(x = gx, y = gy, z = z))
}

## hacked from emdbook

HPDregion <- function(x,vars=1:2,w,h,n=50,prob=0.95,
                          xlab=NULL,ylab=NULL,lims=NULL,...) {
    var1 <- x[,vars[1]]
    var2 <- x[,vars[2]]
    
    varnames <- names(x)[vars]

    if (missing(w)) {
        w <- rep(1/length(var1), length(var1))
    } else {
        if (sum(w) != 1) w <- w/sum(w)
    }
    
    if (is.null(lims)) lims <- c(range(var1), range(var2))
    
    post1 <- kde2d.weighted(var1, var2, w, h, n, lims)
        
    dx <- diff(post1$x[1:2])
    dy <- diff(post1$y[1:2])
    sz <- sort(post1$z)
    ## debugging stuff
    ## if (FALSE) {
    ##   lattice:::contourplot(post1$z)
    ##   d <- with(post1,data.frame(expand.grid(y=y,x=x),z=c(z)))
    ##   lattice:::contourplot(z~x*y,data=d)
    ##   with(post1,contour(x,y,z))
    ##   points(x[,1],x[,2],col=2)
    ##   afun <- function(n) {
    ##     k2 <- kde2d(x[,1],x[,2],n=n,h=c(1,1))
    ##     with(k2,sum(z)*diff(x)[1]*diff(y)[1])
    ##   }
    ##   sapply(5:25,afun)
    ## }
    c1 <- cumsum(sz)*dx*dy
    ## trying to find level containing 95% of volume ...
    levels <- sapply(prob,function(x) {approx(c1,sz,xout=1-x)$y})
    ## meanvec <- c(mean(var1),mean(var2))
    if (is.null(xlab)) xlab <- varnames[1]
    if (is.null(ylab)) ylab <- varnames[2]
    contour(post1$x,post1$y,post1$z,level=levels,
            xlab=xlab,ylab=ylab,drawlabels=FALSE,...)
    invisible(contourLines(post1$x,post1$y,post1$z,levels=levels))
}


## wquant from King et al.
wquant <- function (x, weights, probs = c(0.025, 0.975)) {
    idx <- order(x)
    x <- x[idx]
    weights <- weights[idx]
    w <- cumsum(weights)/sum(weights)
    rval <- approx(w,x,probs,rule=1)
    rval$y
}