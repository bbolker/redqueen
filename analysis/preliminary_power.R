source("../R/powerfun.R")
load("../data/powersim.rda")

samplesize <-seq(20, 200, by=20)

pp <- find_power(simlist, samplesize = samplesize)

image(x=c(5, 10, 15, 20, 25, 30), y=samplesize, pp$power, xlab="Number of sites", ylab="Number of samples per site")
image(x=c(5, 10, 15, 20, 25, 30), y=samplesize, pp$effsize, xlab="Number of sites", ylab="Number of samples per site")

