load("../data/ABC.rda")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf", "psex", "pinf", "psex"
))

param.keep2 <- sapply(res, function(x) sqrt(sum((vergara_CV - x)^2)) < 0.15) & !is.na(param.keep)

print(sum(param.keep2))

par(mfrow=c(1, 3))
ylim <- c(2.5, 6.5, 0.2)
title <- c("meanlog", "sdlog", "bU")

for(i in 1:3) {
    plot(density(param.table[,i]), ylim=c(0, ylim[i]), main=title[i], xlab="")
    lines(density(param.table[param.keep2,i]), col=2)
}

apply(param.table[param.keep2,], 2, mean)
