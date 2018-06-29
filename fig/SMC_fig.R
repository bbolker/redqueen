library(emdbook)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(grid)
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dakr2") scale_fill_brewer(...,palette=palette)

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

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

data_name_simple <- c("D", "M", "V")

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
    sapply(clean_fun, x=comb_smc, simplify=FALSE) %>%
    lapply(mutate, param=1:100)

simlist <- comb_smc %>%
    lapply("[[", "simlist")

SMC_weight <- comb_smc %>%
    lapply("[[", "ww") %>%
    lapply(function(x) lapply(x, function(y) {data.frame(param=1:100, weight=y) })) %>%
    lapply(bind_rows, .id="run") %>%
    bind_rows(.id="fit") %>%
    as.tbl

weighted_list <- clean_list %>%
    lapply(merge, SMC_weight)

SMC_summary <- weighted_list %>%
    lapply(group_by, run, fit, key) %>%
    lapply(summarize, value=weighted.mean(value, w=weight, na.rm=TRUE)) %>%
    lapply(spread, key, value)

if(save) save("comb_summ", "weighted_list", "SMC_weight", "SMC_summary", file="SMC_summary.rda")
if(save) save("comb_summ", file="comb_summ.rda")

gpar <- ggplot(NULL, aes(col=fit, group=run)) +
    geom_line(stat="density", aes(value), lwd=0.9) +
    facet_grid(fit~key, labeller=label_parsed) +
    scale_y_continuous(expand=c(0,0.05)) +
    scale_color_manual(values=c("#D95F02", "#7570B3")) +
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
    filter(fit!="dagan") %>%
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
        labels=data_name[2:3]
    ))

glist <- parlist_factor %>%
    filter(!(key %in% c("asex.genotype"))) %>%
    group_by(key) %>%
    do(plot=gpar%+%.)

names(glist$plot) <- glist$key

glist$plot$`beta[meanlog]` <- glist$plot$`beta[meanlog]` +
    stat_function(fun=function(x) dcauchy(x, location=1, scale=2), col="black", lty=2) +
    xlim(c(-1, 4))

glist$plot$`beta[sdlog]` <- glist$plot$`beta[sdlog]` +
    scale_x_log10(limits=c(0.01, 10)) +
    stat_function(fun=function(x) dlnorm(x, meanlog=0, sdlog=2), col="black", lty=2)

glist$plot$`c[b]` <- glist$plot$`c[b]` +
    xlim(c(0.3,1.8)) +
    stat_function(fun=function(x) dlnorm(x, meanlog=-0.07, sdlog=0.09), col="black", lty=2)

glist$plot$V <- glist$plot$V +
    scale_x_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) +
    stat_function(fun=function(x) dbeta(x, shape1=6, shape2=2), col="black", lty=2)

glist$plot$`epsilon[site]` <- glist$plot$`epsilon[site]`+
    scale_x_continuous(limits=c(0,0.2), breaks=c(0, 0.1)) +
    stat_function(fun=function(x) dbeta(x, shape1=1, shape2=9), col="black", lty=2)

geno_prior <- data.frame(vergara=dbetabinom(0:9, prob=3/9, size=9, theta=5),
           mckone=dbetabinom(0:9, prob=3/9, size=9, theta=5),
           n.genotype=1:10) %>%
    gather(key,value, -n.genotype) %>%
    rename(fit=key) %>%
    mutate(key="G[asex]", fit=factor(fit, labels=data_name[2:3]))

ghist <- parlist_factor %>%
    filter(key=="asex.genotype") %>%
    mutate(key="G[asex]") %>%
    filter(!is.na(value)) %>%
    ggplot() +
        geom_bar(data=geno_prior, aes(n.genotype, value), stat="identity", width=0.5, alpha=0.9) +
        geom_bar(aes(value, y=2*(..count..)/sum(..count..), fill=fit), alpha=0.9, position=position_nudge(x=-0.1), width=0.5) +
        facet_grid(fit~key, labeller=label_parsed) +
        scale_x_continuous(breaks=1:10) +
        scale_y_continuous(expand=c(0,0.05)) +
        scale_fill_manual(values=c("#D95F02", "#7570B3")) +
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
