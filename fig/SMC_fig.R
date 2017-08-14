library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(grid)
library(gridExtra)
load("../data/vergara_summ.rda")
load("../data/dagan_summ.rda")

src <- "../data/"
SMC_files <- c("SMC_vergara_three_loci.rda", "SMC_dagan_three_loci.rda")
summ_files <- c("vergara_summ.rda", "dagan_summ.rda")

get_all <- function(x) {
    ## cat(x,"\n")
    L <- load(x)
    r <- lapply(L,function(x) get(x,envir=parent.frame(2)))
    names(r) <- L
    return(r)
}

comb_smc <- comb_summ <- vector('list', 2)

comb_smc <- lapply(paste0(src, SMC_files), get_all)
comb_summ <- lapply(paste0(src, summ_files), get_all)

names(comb_smc) <- names(comb_summ) <- c("vergara", "dagan")

clean_fun <- function(x, target) {
    x %>%
        lapply("[[", target) %>%
        lapply(bind_rows, .id="run") %>%
        bind_rows(.id="fit") %>%
        as.tbl %>%
        gather(key, value, -run, -fit) %>%
        group_by(run, fit, key)
}  

clean_list <- c("parlist", "sumlist") %>%
    sapply(clean_fun, x=comb_smc, simplify=FALSE)

clean_list %>%
    lapply(summarize, value=mean(value, na.rm=TRUE)) %>%
    lapply(spread, key, value)

gpar <- ggplot(NULL, aes(col=run, group=run)) +
    geom_density(aes(value)) +
    facet_grid(fit~key) +
    theme(
        legend.position = "none",
        strip.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.minor=element_blank()
    )

glist <- clean_list$parlist %>%
    filter(key!="n.genotype") %>%
    group_by(key) %>%
    do(plot=gpar%+%.)

names(glist$plot) <- glist$key

glist$plot$beta.meanlog <- glist$plot$beta.meanlog +
    geom_vline(data=filter(smc_mean, key=="beta.meanlog"), aes(xintercept=mean, col=run), lty=2) +
    stat_function(fun=function(x) dcauchy(x, location=2, scale=1))

glist$plot$beta.sdlog <- glist$plot$beta.sdlog +
    scale_x_log10() +
    stat_function(fun=function(x) dlnorm(x, meanlog=0, sdlog=1))

glist$plot$epsilon.site <- glist$plot$epsilon.site +
    scale_x_log10()

glist$plot$V <- glist$plot$V +
    scale_x_continuous(limits=c(0, 1)) +
    stat_function(fun=function(x) dbeta(x, shape1=6, shape2=2))

ghist <- clean_list$parlist %>%
    filter(key=="n.genotype") %>%
    filter(!is.na(value)) %>%
    group_by(fit, run, key, value) %>%
    summarize(n=n()) %>%
    mutate(freq=n/sum(n)) %>%
    ggplot(aes(value,freq,col=run)) +
        geom_point() +
        geom_line(lty=2) +
        facet_grid(fit~key) +
        theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.minor=element_blank())

gg <- append(glist$plot, list(n.genotype=ghist)) 

do.call(grid.arrange, list(grobs=gg, nrow=1, widths=c(0.8, 0.8, 0.8, 0.8, 1)))

summ_list <- vector('list', 2)

summ_list$dagan <- as.list(dagan_summ)
summ_list$vergara <- as.list(vergara_summ)

summ_df <- summ_list %>%
    lapply(bind_rows) %>%
    lapply(gather) %>%
    bind_rows(.id="fit")

ggplot(clean_list$sumlist) +
    geom_density(aes(value, col=run, group=run)) +
    geom_vline(data=summ_df, aes(xintercept=value)) +
    facet_wrap(fit~key, ncol=6, scale="free")
