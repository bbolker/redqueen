library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())

load("../data/SMC_summary.rda")
load("../data/temporal_effect.rda")
load("../data/true_effect.rda")

temporal_effect <- effect_list %>%
    bind_rows(.id="data") %>%
    gather(key, value, -data, -gen, -sim) %>%
    as.tbl %>%
    group_by(data, sim, key) %>%
    summarize(temporal.effect=mean(value, na.rm=TRUE)) %>%
    rename(test=key)
     
true_effect <- resdf %>%
    filter(test %in% c("spearman", "quantile_quad"), transformation=="raw")  %>%
    mutate(test=factor(test, labels=c("quad", "spearman"))) %>%
    as.tbl %>%
    select(-c(transformation, tgroup)) %>%
    rename(true.effect=effect)

relative_effect <- merge(true_effect, temporal_effect) %>% 
    as.tbl %>%
    group_by %>%
    mutate(relative.effect=temporal.effect/true.effect)

sensitivity <- clean_list$parlist %>%
    filter(run==3) %>%
    mutate(sim=1:n()) %>%
    group_by() %>%
    select(-run) %>%
    rename(param=key, param.value=value, data=fit) %>%
    as.tbl %>%
    merge(relative_effect)
    
ggplot(sensitivity, aes(param.value, relative.effect, col=data)) +
    geom_point() +
    geom_hline(yintercept=1, lty=2) +
    geom_hline(yintercept=0, lty=1) +
    facet_grid(test~param, scale="free") 

ggplot(sensitivity, aes(param.value, temporal.effect, col=data)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(yintercept=0, lty=1) +
    facet_grid(test~param, scale="free") 

ggplot(sensitivity, aes(param.value, true.effect, col=data)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_hline(yintercept=0, lty=1) +
    facet_grid(test~param, scale="free") 

ggplot(sensitivity, aes(true.effect, temporal.effect)) +
    geom_point() +
    geom_hline(yintercept=0, lty=2) +
    geom_vline(xintercept=0, lty=2) +
    geom_abline(lty=1) +
    facet_wrap(test~data, scale="free", nrow=2)

