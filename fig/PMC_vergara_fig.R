library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(grid)
library(gridExtra)
load("../data/PMC_vergara_three_loci.rda")
load("../data/vergara_summ.rda")

pardf <- parlist %>%
    bind_rows(.id="tol") %>%
    mutate(tol=factor(tol, labels = c(0.8, 0.4, 0.2))) %>%
    mutate(beta.sdlog=exp(log.beta.sdlog),
           V=plogis(logit.V)) %>%
    select(tol, beta.meanlog, beta.sdlog, V) %>%
    gather(key, value, -tol) %>%
    group_by(tol, key)

parmean <- pardf %>%
    summarize(value=mean(value, na.rm=TRUE))

gpar <-ggplot(pardf) +
    geom_density(aes(value, col=tol)) +
    facet_wrap(~key, scale="free")

gpar +
    geom_vline(data=parmean, aes(xintercept=value, col=tol), lty=2)

sumdf <- sumlist %>%
    bind_rows(.id="tol") %>%
    mutate(tol=factor(tol, labels = c(0.8, 0.4, 0.2))) %>%
    gather(key, value, -tol) %>%
    group_by(tol, key)

vergdf <- data.frame(
    key=names(vergara_summ),
    value=vergara_summ
)

gpar %+% sumdf +
    geom_vline(data=vergdf,aes(xintercept=value), lty=2)


