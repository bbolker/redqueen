library(ggplot2)
library(tidyr)
library(dplyr)
source("../R/ABC_funs.R")
load("../data/vergara_summ.rda")
load("../data/dagan_summ.rda")
load("../data/fitted_sim.rda")

simlist <- list(
    dagan=dagan_sim,
    vergara=vergara_sim
)

sites <- c(22, 4)

summ_list <- vector('list', 2)

for(i in 1:2) {
    siml <- simlist[[i]]
    rl <- vector('list', length(siml))
    for(j in 1:length(siml)) {
        sim <- siml[[j]]
        rl[[j]] <- as.data.frame(t(replicate(10, sumfun(sim, sitesample=sites[i]))))
    }
    summ_list[[i]] <- do.call('rbind', rl)
}

names(summ_list) <- names(simlist)

summ_df <- summ_list %>%
    bind_rows(.id="data") %>%
    gather(key, value, -data)

mean_df <- summ_df %>%
    group_by(data, key) %>%
    summarize(
        mean=mean(value)
    )

orig_df <- list(dagan=dagan_summ,
     vergara=vergara_summ) %>%
    lapply(as.list) %>%
    lapply(as.data.frame) %>%
    bind_rows(.id="data") %>%
    gather(key, value, -data)

ggplot(summ_df) +
    geom_density(aes(value)) +
    geom_vline(data=orig_df, aes(xintercept=value), lty=2, col="red") +
    geom_vline(data=mean_df, aes(xintercept=mean), lty=2) +
    facet_wrap(data~key, scale="free", ncol=6)

