load("../data/ABC.rda")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf", "psex", "pinf", "psex"
))

param.keep <- sapply(res, function(x) sqrt(sum((vergara_CV - x)^2)) < 0.18)

par(mfrow=c(1, 3))
ylim <- c(2, 5, 0.2)

for(i in 1:3) {
    plot(density(param.table[,i]), ylim=c(0, ylim[i]))
    lines(density(param.table[param.keep,i]), col=2)
}

apply(param.table[param.keep,], 2, mean)
