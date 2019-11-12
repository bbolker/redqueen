library(ggplot2); theme_set(theme_bw())
library(scales) ## for trans_new?
library(dplyr)
library(tidyr)
library(readxl)
library(quantreg)

dagan <- read_excel("../data/JEB12245-Dryad.xlsx", sheet=1)

dagan_summ <- dagan %>%
    select_("Males", "Infected") %>%
    mutate(sex=2*Males) %>%
    summarize_all(funs(mean, sd(.)/mean(.))) %>%
    select(5, 6, 2, 3) %>%
    setNames(c("pinf.siteCV", "psex.siteCV", "pinf.mean", "psex.mean")) %>%
    unlist

save("dagan_summ", file="dagan_summ.rda")
 
logitfun <- function(x, min.p) log((x+min.p)/(1-x+min.p))

dagan_trans <- dagan %>%
    mutate(
        logit.m=logitfun(Males, 1/100),
        logit.i=logitfun(Infected, 2/35)
    ) %>%
    select(PopulationName, logit.m, logit.i)

lm1 <- lm(logit.m~logit.i, data=dagan_trans)

summary(lm1); AIC(lm1)

lm2 <- lm(logit.m~poly(logit.i, degree=2), data=dagan_trans)

summary(lm2); AIC(lm2)

rr1 <- rq(logit.m~logit.i, data=dagan_trans, tau=seq(0.1, 0.9, 0.1))

plot(rr1)

summary(rr1)[[5]]; AIC(rr1)[[5]]

rr2 <- rq(logit.m~poly(logit.i, degree=2), data=dagan_trans, tau=0.5)

summary(rr2); AIC(rr2)

## par(mfrow=c(1,1))
## plot(dagan_trans[,-1])

newlogit <- seq(min(dagan_trans$logit.i), max(dagan_trans$logit.i), by=0.1)

prr <- predict(rr2, 
               newdata=data.frame(logit.i=newlogit), 
               interval="confidence") %>%
    as.data.frame %>%
    mutate(logit.i=newlogit) %>%
    rename(logit.m=fit)

inverse_logit_trans <- trans_new("inverse logit",
    transform = plogis,
    inverse = qlogis)

ggplot(dagan_trans, aes(logit.i, logit.m)) +
    geom_point() +
    geom_quantile(quantiles=0.5, formula=y~poly(x, degree=2), col=2, lwd=1) +
    geom_smooth(method='lm', formula=y~x, col="blue", fill="blue", alpha=0.1, lwd=1) +
    geom_ribbon(data=prr, aes(logit.i, ymin=lower, ymax=higher), alpha=0.1, fill=2)

## lin reg

ggplot(dagan_trans, aes(logit.i, logit.m)) +
    geom_point() +
    geom_quantile(quantiles=seq(0.05, 0.95, 0.1), formula=y~poly(x, degree=2), col=2, lwd=1) +
    geom_smooth(method='lm', formula=y~poly(x, degree=2), col="blue", fill="blue", alpha=0.1, lwd=1)
