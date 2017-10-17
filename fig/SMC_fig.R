library(emdbook)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_size = 14,
                                     base_family = "Times"))
library(grid)
library(gridExtra)

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

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

data_name <- c(expression(Dagan~italic(et~al.)~"(2005)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

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
    geom_line(stat="density", aes(value)) +
    facet_grid(fit~key, labeller=label_parsed) +
    scale_y_continuous(expand=c(0,0.05)) +
    theme(
        legend.position = "none",
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        strip.text.y=element_blank(),
        axis.title=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        plot.margin=unit(c(1,0.1,1,0.1), "cm")
    )

parlist_factor <- clean_list$parlist %>%
    group_by %>%
    filter(run==3) %>%
    mutate(key=factor(key,
        labels=c(
            "beta[meanlog]",
            "beta[sdlog]",
            "c[b]",
            "epsilon[site]",
            "asex.genotype",
            "V"
        )                  
    ), fit=factor(fit,
        labels=data_name
    ))

glist <- parlist_factor %>%
    filter(!(key %in% c("asex.genotype"))) %>%
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
    mutate(key="beta[meanlog]", run="1", fit=factor(fit, labels=data_name))

glist$plot$`beta[meanlog]` <- glist$plot$`beta[meanlog]` +
    geom_line(data=beta_prior, aes(x, value), col=1) +
    xlim(c(-2, 4))

glist$plot$`beta[sdlog]` <- glist$plot$`beta[sdlog]` +
    scale_x_log10(limits=c(1e-3, 10)) +
    stat_function(fun=function(x) dlnorm(x, meanlog=0, sdlog=2), col="black")

glist$plot$`c[b]` <- glist$plot$`c[b]` +
    xlim(c(0.3,1.8)) +
    stat_function(fun=function(x) dlnorm(x, meanlog=-0.1, sdlog=0.1), col="black")

glist$plot$V <- glist$plot$V +
    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
    stat_function(fun=function(x) dbeta(x, shape1=6, shape2=2), col="black")

epsilon_x <- exp(seq(log(1e-5), 1, 0.01))

epsilon_prior <- data.frame(vergara=dbeta(epsilon_x, shape1=1, shape2=9),
                         dagan=dbeta(epsilon_x, shape1=1, shape2=9),
                         mckone=dbeta(epsilon_x, shape1=2, shape2=8),
                         x=epsilon_x) %>%
    gather(key,value, -x) %>%
    rename(fit=key) %>%
    mutate(key="epsilon[site]", run="1", fit=factor(fit, labels=data_name))

glist$plot$`epsilon[site]` <- glist$plot$`epsilon[site]`+
    geom_line(data=epsilon_prior, aes(x, value), col=1) +
    scale_x_continuous(limits=c(0,0.2))

geno_prior <- data.frame(vergara=dbetabinom(0:9, prob=2/9, size=9, theta=5),
           dagan=dbetabinom(0:9, prob=5/9, size=9, theta=5),
           mckone=dbetabinom(0:9, prob=2/9, size=9, theta=5),
           n.genotype=1:10) %>%
    gather(key,value, -n.genotype) %>%
    rename(fit=key) %>%
    mutate(key="G[asex]", fit=factor(fit, labels=data_name))

ghist <- parlist_factor %>%
    filter(key=="asex.genotype") %>%
    mutate(key="G[asex]") %>%
    filter(!is.na(value)) %>%
    ggplot() +
        geom_bar(data=geno_prior, aes(n.genotype, value), stat="identity", width=0.5) +
        geom_bar(aes(value, y=3*(..count..)/sum(..count..), fill=run), position=position_nudge(x=-0.1), width=0.5) +
        facet_grid(fit~key, labeller=label_parsed) +
        scale_x_continuous(breaks=1:10) +
        scale_y_continuous(expand=c(0,0.05)) +
        theme(axis.title=element_blank(),
            legend.position = "none",
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"),
            panel.grid=element_blank(),
            panel.spacing=grid::unit(0,"lines"),
            plot.margin=unit(c(1,0.3,1,0.1), "cm"))

glist$plot <- lapply(glist$plot, function(x) {x$layers <- rev(x$layers); x})

gg <- append(glist$plot, list(`asex[genotype]`=ghist))

gg_smc_param <- do.call(arrangeGrob, list(grobs=gg, nrow=1, widths=c(0.9, 0.9, 0.9, 0.9, 0.9, 1.1)))

if (save) ggsave("smc_param.pdf", gg_smc_param, width=8, height=5)

gg_reduced <- gg[c("beta[meanlog]", "beta[sdlog]", "c[b]", "V")] %>%
    lapply(function(g) g +  theme(plot.margin=unit(c(0.2,0.05,0.2,0.05), "cm")))

gg_reduced$V <- gg_reduced$V +
    theme(plot.margin=unit(c(0.2,0.5,0.2,0.05),"cm"),
          strip.text.y=element_text(angle=-90, hjust=1))

gg_smc_param_reduced <- do.call(arrangeGrob, list(grobs=gg_reduced, nrow=1, widths=c(0.9, 0.9, 0.9, 1)))

if (save) ggsave("smc_param_red.pdf", gg_smc_param_reduced, width=8, height=5)

summ_df <- comb_summ  %>%
    lapply(function(x) as.list(x[[1]])) %>%
    lapply(bind_rows) %>%
    lapply(gather) %>%
    bind_rows(.id="fit") %>%
    mutate(fit=factor(fit, labels=data_name)) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean", "across~site~CV", "across~generation~CV"),2)))

gg_summary_df <- clean_list$sumlist %>% 
    filter(run==3) %>%
    group_by() %>%
    mutate(fit=factor(fit, labels=data_name)) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean", "across~site~CV", "across~generation~CV"),2)))

gg_summary_mean <- SMC_summary$sumlist %>% 
    group_by() %>%
    mutate(fit=factor(fit, labels=data_name)) %>%
    filter(run==3) %>%
    gather(key, value, -fit, -run) %>%
    mutate(gvar=ifelse(grepl("pinf", key), "proportion~infected", "proportion~sexual")) %>%
    mutate(key=factor(key,
                      levels=c("pinf.mean", "pinf.siteCV", "pinf.timeCV", "psex.mean", "psex.siteCV", "psex.timeCV"),
                      labels=rep(c("mean", "across~site~CV", "across~generation~CV"),2)))

blank_df <- data.frame(
    fit=gg_summary_mean$fit[1],
    key="mean",
    value=0,
    gvar="proportion~infected"
)

gg_smc_summ <- ggplot(gg_summary_df) +
    geom_blank(data=blank_df, aes(value)) +
    geom_line(stat="density", aes(value, col=run, group=run)) +
    geom_vline(data=gg_summary_mean, aes(xintercept=value, col=run, group=run)) +
    geom_vline(data=summ_df, aes(xintercept=value)) +
    facet_grid(fit~gvar+key, scale="free", labeller = label_parsed) +
    scale_x_continuous(name="") +
    theme(
        legend.position = "none"
    )

if (save) ggsave("smc_summary.pdf", gg_smc_summ, width=8, height=6)

