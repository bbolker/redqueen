library(ggplot2); theme_set(theme_bw())
library(dplyr)
load("../data/fitted_power.rda")

level <- 0.05

simdf <- simlist %>%
    lapply(bind_rows, .id="test") %>%
    bind_rows(.id="sim") %>%
    mutate(sites=as.factor(sites)) %>%
    as.tbl

sumdf <- simdf %>%
    group_by(sim, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

sig_sumdf <- simdf %>%
    filter(p.value < level) %>%
    group_by(sim, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

sig_simdf <- simdf %>%
    group_by(sim, test, sites, samples) %>%
    summarize(power=mean(p.value<=level, na.rm=TRUE))
    
(gg_sum <- ggplot(sumdf, aes(samples, median, fill=sites)) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), alpha=0.3) +
    geom_point(aes(col=sites)) +
    geom_line(aes(col=sites)) +
    facet_grid(test~sim, scale="free") +
    geom_hline(yintercept=0, lty=2))

gg_sum %+% sig_sumdf

ggplot(sig_simdf, aes(samples, power, col=sites)) +
    geom_point() +
    geom_line(lty=2) +
    facet_grid(sim~test)

