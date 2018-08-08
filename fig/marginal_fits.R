library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
source("../R/HPDregion.R")
source("../R/ABC_funs.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

source("../R/HPDregion.R")

load_names <- c("dagan", "mckone", "vergara")

dir <- "../simdata/"
slist <- vector('list', 3)
names(slist) <- load_names

for (n in load_names) {
    ll <- vector('list', 10)
    
    for (i in 0:9) {
        fn <- paste0("post_", n, "_summ_",i, ".rda")
        load(paste0(dir, fn))
        ll[[i + 1]] <- do.call("rbind", summ_reslist)
    }
    
    slist[[n]] <- do.call("rbind", ll)
}

load("../data/SMC_summary.rda")

accepted_df <- weighted_list$sumlist

load("../data/comb_summ.rda")

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

names(slist) <- data_name

gg_summary_df <- slist %>%
    bind_rows(.id="fit") %>%
    as.tbl %>%
    gather(key, value, -sim, -param, -fit, -weight) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean~proportion", "CV~across~population", "CV~across~generation"),2))) %>%
    filter(!(fit!=as.character(data_name[3]) & key=="CV~across~generation")) %>%
    mutate(weight=weight/1000)

gg_summary_quant <- gg_summary_df %>%
    group_by(fit, key, gvar) %>%
    summarise(
        lwr=wquant(value, weight, 0.025),
        upr=wquant(value, weight, 0.975)
    )

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

gg_accepted <- accepted_df %>%
    filter(run==4, !is.na(value)) %>%
    mutate(fit=factor(fit, labels=data_name)) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean~proportion", "CV~across~population", "CV~across~generation"),2)))

gg_smc_summ <- ggplot(gg_summary_quant, aes(fit)) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1) +
    geom_violin(data=gg_accepted, aes(y=value, fill=fit, weight=weight), alpha=0.5, width=0.5) +
    geom_errorbar(data=summ_df, aes(ymin=value, ymax=value), lty=2) +
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

gg_summary_df %>%
    group_by(fit, key, gvar) %>%
    filter(key %in% c("mean~proportion")) %>%
    summarize(mean=weighted.mean(value, weight=weight),
              lwr=wquant(value, weight=weight, 0.025),
              upr=wquant(value, weight=weight, 0.975))
