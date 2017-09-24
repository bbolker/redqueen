library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(dplyr)
load("../data/fitted_power.rda")
load("../data/SMC_summary.rda")

level <- 0.05

reslist2 <- reslist # %>%
#    lapply(function(x) mut2ate(x, p.value=ifelse(is.na(p.value), 1, p.value),
#                              effect.size=ifelse(is.na(effect.size), 0, effect.size))) %>%
#    lapply(filter, !grepl("Error", effect.size)) ## remove try-error    

simdf <- reslist2 %>%
    bind_rows(.id="data") %>%
    as.tbl %>%
    mutate(sites=factor(sites))
    
sig_simdf <- simdf %>%
    group_by(data, test, sites, samples, sim) %>%
    summarize(positive.power=mean(p.value<=level & effect.size > 0, na.rm=TRUE),
              negative.power=mean(p.value<=level & effect.size < 0, na.rm=TRUE)) %>%
    gather(key, value, -data, -test, -sites, -samples, -sim) %>%
    group_by(data, test, sites, samples, key) %>%
    summarize(median=median(value),
              lwr=quantile(value, 0.025, na.rm=TRUE),
              upr=quantile(value, 0.975, na.rm=TRUE)) %>%
    filter(sites %in% c(10, 20, 30))

ggplot(sig_simdf, aes(samples, median, group=sites, col=sites)) +
    geom_point() +
    geom_ribbon(aes(ymin=lwr, ymax=upr, fill=sites), alpha=0.1) +
    geom_line() +
    facet_grid(data~test+key)

sumdf <- simdf %>%
    group_by(data, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

sig_sumdf <- simdf %>%
    filter(p.value < level) %>%
    group_by(data, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

(gg_sum <- ggplot(sumdf, aes(samples, median, fill=sites)) +
        geom_ribbon(aes(ymax=upr, ymin=lwr), alpha=0.3) +
        geom_point(aes(col=sites)) +
        geom_line(aes(col=sites)) +
        facet_grid(test~data, scale="free") +
        geom_hline(yintercept=0, lty=2))

gg_sum %+% sig_sumdf

pardf <- clean_list$parlist %>%
    filter(run==3) %>%
    mutate(sim=1:50) %>%
    rename(data=fit)

alldf <- simdf %>% 
    group_by(data,sim,sites,samples, test) %>%
    summarize(
        mean.effect=mean(effect.size, na.rm=TRUE),
        positive.power=mean(p.value<level & effect.size > 0, na.rm=TRUE),
        negative.power=mean(p.value<level & effect.size < 0, na.rm=TRUE)
    ) %>%
    gather(key, value, -data, -sim, -sites, -samples, -test, -mean.effect) %>%
    rename(type=key, power=value) %>%
    merge(pardf) %>%
    filter(samples==100, sites==20)

ggplot(alldf, aes(value, mean.effect, col=data)) +
    geom_point(alpha=0.3) +
    geom_smooth(method='lm') +
    facet_grid(test~key, scale="free") +
    geom_hline(yintercept=0, lty=2)

pardf_mean <- pardf %>%
    summarize(mean=mean(value))

ggplot(alldf, aes(value, power, col=data)) +
    geom_point(alpha=0.3) +
    geom_smooth(method='lm', aes(fill=data), alpha=0.2) +
    geom_vline(data=pardf_mean, aes(xintercept=mean, col=data), lty=2) + 
    scale_y_continuous(limits=c(0, 1)) +
    facet_grid(test+type~key, scale="free_x")
