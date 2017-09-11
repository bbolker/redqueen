library(emdbook)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(grid)
library(gridExtra)

save <- FALSE

src <- "../data/"
SMC_files <- c("SMC_vergara.rda", "SMC_dagan.rda", "SMC_mckone.rda")
summ_files <- c("vergara_summ.rda", "dagan_summ.rda", "mckone_summ.rda")

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

names(comb_smc) <- names(comb_summ) <- c("vergara", "dagan", "mckone")

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

simlist <- comb_smc %>%
    lapply("[[", "simlist")

SMC_summary <- clean_list %>%
    lapply(summarize, value=mean(value, na.rm=TRUE)) %>%
    lapply(spread, key, value)

if(save) save("comb_smc", "clean_list", "simlist", "SMC_summary", file="SMC_summary.rda")

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

beta_x <- seq(-3, 5, 0.01)

beta_prior <- data.frame(vergara=dcauchy(beta_x, location=2, scale=1),
                         dagan=dcauchy(beta_x, location=0, scale=2),
                         mckone=dcauchy(beta_x, location=0, scale=2),
                         x=beta_x) %>%
    gather(key,value, -x) %>%
    rename(fit=key) %>%
    mutate(key="beta.meanlog", run=1)

glist$plot$beta.meanlog <- glist$plot$beta.meanlog +
    geom_line(data=beta_prior, aes(x, value), col=1) +
    xlim(c(-2, 4))

glist$plot$beta.sdlog <- glist$plot$beta.sdlog +
    scale_x_log10() +
    stat_function(fun=function(x) dlnorm(x, meanlog=0, sdlog=1), col="black")

epsilon_x <- exp(seq(log(1e-10), log(3e-1), 0.01))

epsilon_prior <- data.frame(vergara=dbeta(epsilon_x, shape1=1, shape2=99),
                         dagan=dbeta(epsilon_x, shape1=1, shape2=99),
                         mckone=dbeta(epsilon_x, shape1=1, shape2=19),
                         x=epsilon_x) %>%
    gather(key,value, -x) %>%
    rename(fit=key) %>%
    mutate(key="epsilon.site", run=1)

glist$plot$epsilon.site <- glist$plot$epsilon.site +
    geom_line(data=epsilon_prior, aes(x, value), col=1)

glist$plot$c_b <- glist$plot$c_b +
    xlim(c(0.3,1.8)) +
    stat_function(fun=function(x) dlnorm(x, meanlog=-0.1, sdlog=0.1), col="black")

glist$plot$V <- glist$plot$V +
    scale_x_continuous(limits=c(0, 1)) +
    stat_function(fun=function(x) dbeta(x, shape1=6, shape2=2), col="black")

geno_prior <- data.frame(vergara=dbetabinom(0:9, prob=2/9, size=9, theta=5),
           dagan=dbetabinom(0:9, prob=5/9, size=9, theta=5),
           mckone=dbetabinom(0:9, prob=2/9, size=9, theta=5),
           n.genotype=1:10) %>%
    gather(key,value, -n.genotype) %>%
    rename(fit=key) %>%
    mutate(key="n.genotype")

ghist <- clean_list$parlist %>%
    filter(key=="n.genotype") %>%
    filter(!is.na(value)) %>%
    ggplot() +
        geom_histogram(aes(value, ..density.., fill=run), position="identity", alpha=0.5, bins=10) +
        geom_line(data=geno_prior, aes(n.genotype, value)) +
        facet_grid(fit~key) +
        theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.minor=element_blank())

gg <- append(glist$plot, list(n.genotype=ghist))

do.call(grid.arrange, list(grobs=gg, nrow=1, widths=c(0.8, 0.8, 0.8, 0.8, 0.8, 1)))

summ_df <- comb_summ  %>%
    lapply(function(x) as.list(x[[1]])) %>%
    lapply(bind_rows) %>%
    lapply(gather) %>%
    bind_rows(.id="fit")

ggplot(clean_list$sumlist) +
    geom_density(aes(value, col=run, group=run)) +
    geom_vline(data=summ_df, aes(xintercept=value)) +
    geom_vline(data=gather(SMC_summary$sumlist, key, value, -fit, -run), aes(xintercept=value, col=run, group=run), lty=2) +
    facet_wrap(fit~key, ncol=6, scale="free")
