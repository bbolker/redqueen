library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

source("../R/ABC_funs.R")

load("../data/SMC_summary.rda")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

nsim <- 100

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

simlist <- list(
    dagan=simlist$dagan[[4]],
    mckone=simlist$mckone[[4]],
    vergara=simlist$vergara[[4]]
)

sites <- c(22, 18, 4)

reslist <- vector('list', length(simlist))

names(simlist) <- names(reslist) <- data_name

set.seed(101)
for(i in 1:length(simlist)) {
    sim <- simlist[[i]]
    
    nvalid <- sum(sapply(sim, length)>0)
    
    subreslist <- vector('list', nvalid)
    
    for(j in 1:nvalid) {
        s <- sim[[j]]
        subreslist[[j]] <- as.data.frame(t(replicate(nsim, sumfun(s, sitesample=sites[i]))))
    }
    reslist[[i]] <- subreslist
}

gg_summary_df <- reslist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="fit") %>%
    as.tbl %>%
    gather(key, value, -sim, -fit) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean~proportion", "CV~across~population", "CV~across~generation"),2))) %>%
    filter(!(fit!=as.character(data_name[3]) & key=="CV~across~generation"))

summ_df <- comb_summ  %>%
    lapply(function(x) as.list(x[[1]])) %>%
    lapply(bind_rows) %>%
    lapply(gather) %>%
    bind_rows(.id="fit") %>%
    mutate(fit=factor(fit, labels=data_name)) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean~proportion", "CV~across~population", "CV~across~generation"),2)))

gg_smc_summ <- ggplot(gg_summary_df, aes(fit, value)) +
    geom_violin(aes(fill=fit), alpha=0.5, width=0.7) +
    geom_point(data=summ_df, shape=0, size=2.5) +
    facet_grid(key~gvar, scale="free", labeller = label_parsed) +
    scale_fill_discrete(label=data_name) +
    theme(
        strip.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x=element_blank(),
        axis.title=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position="bottom",
        legend.title = element_blank()
    )

if (save) ggsave("smc_summary.pdf", gg_smc_summ, width=6, height=6)

reslist %>%
    lapply(bind_rows) %>%
    bind_rows(.id="data") %>%
    group_by(data) %>%
    summarise_all(mean, na.rm=TRUE)

comb_summ %>%
    lapply(as.data.frame) %>%
    lapply(t) %>%
    lapply(as.data.frame) %>%
    bind_rows(.id="data")
