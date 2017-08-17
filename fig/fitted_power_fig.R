library(ggplot2)
library(dplyr)
load("../data/fitted_power.rda")

level <- 0.05

simdf <- simlist %>%
    lapply(function(x) do.call("rbind", x)) %>%
    lapply(function(x) apply(x, 1, function(x) bind_rows(x, .id="sites"))) %>%
    lapply(bind_rows, .id="test") %>%
    bind_rows(.id="sim") %>%
    as.tbl %>%
    mutate(sites=factor(sites, label=seq(10, 50, by=5)))

sig_simdf <- simdf %>%
    group_by(sim, test, sites) %>%
    summarize(power=mean(p.value<=level, na.rm=TRUE))
    
ggplot(simdf, aes(sites, effect.size)) +
    geom_violin() +
    facet_grid(test~sim, scale="free") +
    geom_hline(yintercept=0, lty=2)

ggplot(sig_simdf) +
    geom_point(aes(sites, power)) +
    geom_line(aes(sites, power), lty=2) +
    facet_grid(sim~test)
