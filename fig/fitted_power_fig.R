library(ggplot2); theme_set(theme_bw())
library(tidyr)
library(dplyr)
load("../data/fitted_power.rda")
load("../data/SMC_summary.rda")

level <- 0.05

reslist2 <- reslist %>%
    lapply(bind_rows, .id="test") %>%
    lapply(filter, !grepl("Error", effect.size)) ## remove try-error    

simdf <- lapply(reslist2, function(x){
    df <- as.data.frame(sapply(x[,-1], as.numeric))
    df$test <- x[,1]
    df
}) %>%
    bind_rows(.id="data") %>%
    as.tbl %>%
    mutate(sites=factor(sites))
    
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

sig_simdf <- simdf %>%
    group_by(data, test, sites, samples) %>%
    summarize(power=mean(p.value<=level, na.rm=TRUE))

(gg_sum <- ggplot(sumdf, aes(samples, median, fill=sites)) +
    geom_ribbon(aes(ymax=upr, ymin=lwr), alpha=0.3) +
    geom_point(aes(col=sites)) +
    geom_line(aes(col=sites)) +
    facet_grid(test~data, scale="free") +
    geom_hline(yintercept=0, lty=2))

gg_sum %+% sig_sumdf

ggplot(sig_simdf, aes(samples, power, col=sites)) +
    geom_point() +
    geom_line(lty=2) +
    facet_grid(data~test)

pardf <- clean_list$parlist %>%
    filter(run==2) %>%
    mutate(sim=1:50) %>%
    rename(data=fit)

alldf <- simdf %>% 
    group_by(data,sim,sites,samples, test) %>%
    summarize(
        median.effect=median(effect.size, na.rm=TRUE),
        power=mean(p.value<level, na.rm=TRUE)
    ) %>%
    merge(pardf) %>%
    filter(samples==100, sites==20)

ggplot(alldf, aes(value, median.effect, col=data)) +
    geom_point() +
    geom_smooth(method='lm') +
    facet_grid(test~key, scale="free") +
    geom_hline(yintercept=0, lty=2)

ggplot(alldf, aes(value, power, col=data)) +
    geom_point() +
    geom_smooth(method='lm') +
    scale_y_continuous(limits=c(0, 1)) +
    facet_grid(test~key, scale="free_x")
