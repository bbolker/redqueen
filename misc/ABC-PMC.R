n <- 2000
tmax <- 5

samp <- rnorm(n)

weight <- rep(1/n, n)

theta.list <- vector('list', tmax)
theta.list[[1]] <- samp

weight.list <- vector('list', tmax)

V <- 2 * var(samp)

for (i in 2:tmax) {
    samp <- sample(theta.list[[i-1]], prob=weight, replace=TRUE)
    
    theta <- rnorm(n, mean=samp, sd=sqrt(V))
    
    theta.list[[i]] <- theta
    
    weight <- dnorm(theta)/(
        sapply(theta, function(x) sum(dnorm(theta.list[[i-1]], mean=x, sd=sqrt(V)) * weight))
    )
    
    weight.list[[i]] <- weight <- weight/sum(weight)
    
    V <- 2 * cov.wt(as.data.frame(theta), wt=weight)$cov
}

plot(density(theta.list[[1]]))
lines(density(theta.list[[2]], weight=weight.list[[2]]), col=2)
lines(density(theta.list[[3]], weight=weight.list[[3]]), col=3)
lines(density(theta.list[[5]], weight=weight.list[[5]]), col=4)

plot(density(theta.list[[5]], weight=weight.list[[5]]))
curve(dnorm(x), add=TRUE, col=2)
