library(dplyr)
library(ggplot2)
source("../R/powerfun.R")
load("../data/beta_sim.rda")

level <- 0.05

tests <- list(
    linear=test_lm,
    quadratic=test_quad
)

power_list <- list()

for (test_name in names(tests)) {
    res <- lapply(simlist, powerfun, nsim=10, nsample=100, test=tests[[test_name]])
    power_list[[test_name]] <- bind_rows(res, .id="beta.meanlog")
}

effect_df <- power_list %>%
    bind_rows(.id="test") %>%
    as.tbl

sig_effect_df <- effect_df %>%
    filter(p.value < level)

power_df <- effect_df %>%
    group_by(test, beta.meanlog) %>%
    summarize(power=mean(p.value<level))

ggplot(effect_df, aes(beta.meanlog, effect.size)) +
    geom_violin(fill=NA) +
    geom_violin(data=sig_effect_df, col=2, fill=NA) +
    geom_hline(yintercept=0, lty=2) +
    facet_grid(test~.)
    
ggplot(power_df, aes(as.numeric(beta.meanlog), power)) +
    geom_point() +
    geom_line() +
    facet_grid(test~.)
