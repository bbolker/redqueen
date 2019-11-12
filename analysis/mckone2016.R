library(ggplot2); theme_set(theme_bw())
library(scales) ## for trans_new?
library(dplyr)
library(tidyr)
library(quantreg)

mckone <- read.csv("../data/mckone2014.csv", header=FALSE)

mckone <- (mckone/100) %>%
    setNames(c("infection", "male")) %>%
    mutate(male=round(male, 3),
           infection=round(infection, 3))

mckone_summ <- mckone %>%
    mutate(sex=2*male) %>%
    summarize_all(funs(mean, sd(.)/mean(.))) %>%
    select(4, 6, 1, 3) %>%
    setNames(c("pinf.siteCV", "psex.siteCV", "pinf.mean", "psex.mean")) %>%
    unlist

save("mckone_summ", file="mckone_summ.rda")

logitfun <- function(x, min.p) log((x+min.p)/(1-x+min.p))

mckone_trans <- mckone %>%
    mutate(
        logit.m=logitfun(male, 1/100),
        logit.i=logitfun(infection, 1/100)
    ) %>%
    select(logit.m, logit.i)

lm1 <- lm(logit.m~logit.i, data=mckone_trans)

summary(lm1); AIC(lm1)

lm2 <- lm(logit.m~poly(logit.i, degree=2), data=mckone_trans)

summary(lm2); AIC(lm2)

ggplot(mckone_trans, aes(logit.i, logit.m)) +
    geom_point() +
    geom_smooth(method="lm")


