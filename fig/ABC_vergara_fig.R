library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

load("../data/ABC_vergara.rda")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf.timeCV", "psex.timeCV", "pinf.siteCV", "psex.siteCV"
    ))

vergara_mean <- structure(
    c(0.433838731993719, 0.703887586155451), 
    .Names = c("pinf.mean", 
               "psex.mean"))

vergara_summ <- c(vergara_CV, vergara_mean)

new.tol <- 0.36 ## sqrt(6*0.15^2)

## sum(sapply(res, inherits, "try-error"))

param.keep2 <- sapply(res, function(x) ifelse(is.null(x) || inherits(x, "try-error"), FALSE, sqrt(sum((vergara_summ - x)^2)) < new.tol))

sum(param.keep2)

param.table2 <- cbind(param.keep2, param.table)

prior <- param.table2 %>%
    as_data_frame %>%
    gather(key, value, -param.keep2)

posterior <- prior %>%
    filter(param.keep2==TRUE) %>%
    select(-param.keep2)

ggplot(prior, aes(x=value)) +
    geom_density() +
    geom_density(data=posterior, col=2) +
    facet_wrap(~key, scales="free")

posterior %>%
    group_by(key) %>%
    summarise(mean=mean(value))


