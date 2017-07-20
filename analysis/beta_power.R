source("../R/powerfun.R")
load("../data/betasim.rda")

meanlog <- seq(0.2, 3, by=0.4)
samplesize <-seq(20, 200, by=20)

bb <- find_power(simlist, samplesize = samplesize)

contour(x=meanlog, y=samplesize, bb$power, xlab="mean log beta", ylab="Number of samples per site")
contour(x=meanlog, y=samplesize, bb$effsize, xlab="mean log beta", ylab="Number of samples per site")

ss <- with(as.list(simlist[3,5][[1]]), {
    S <- S.count
    A <- A.count
    I <- SI.count + AI.count
    N <- S+A
    list(S=S/N, A=A/N, I=I/N, SI=SI.count/N, AI=AI.count/N)
})

plot(ss$I[1000,], ss$S[1000,])

tsub <- seq(1000, 1100, by=5)
matplot((1-ss$SI)/(1-ss$AI), type="l", xlim=c(min(tsub), max(tsub)), ylim=c(0, 10))
plot((1-ss$SI[tsub,3])/(1-ss$AI[tsub,3]), type="l")

aa <- acf((1-ss$SI[tsub,3])/(1-ss$AI[tsub,3]))

plot((1-rowMeans(ss$SI[tsub,]))/(1-rowMeans(ss$AI[tsub,])), type="l")
acf((1-rowMeans(ss$SI[tsub,]))/(1-rowMeans(ss$AI[tsub,])))



