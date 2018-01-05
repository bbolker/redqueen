library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)
load("../data/fitted_power.rda")
load("../data/SMC_summary.rda")

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
}

save <- FALSE

data_name <- c(expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

level <- 0.05

simdf <- reslist %>%
    bind_rows(.id="data") %>%
    as.tbl %>%
    mutate(sites=factor(sites),
           p.value=ifelse(is.na(p.value), 1, p.value),
           effect.size=ifelse(is.na(effect.size), 0, effect.size))

spread_simdf <- simdf %>%
    group_by(sim, test) %>%
    mutate(id=1:n()) %>%
    gather(key, value, -data, -sites, -test, -sim, -samples, -id) %>%
    spread(test, value)

sig_simdf <- simdf %>%
    group_by(data, test, sites, samples) %>%
    summarize(positive.power=mean(p.value<=level & effect.size > 0),
              negative.power=mean(p.value<=level & effect.size < 0)) %>%
    gather(key, value, -data, -test, -sites, -samples) %>%
    group_by() %>%
    mutate(key=factor(key, labels=c("negative~correlation", "positive~correlation")))

spearman_power <- sig_simdf %>%
    filter(test=="spearman") %>%
    mutate(data=factor(data, labels=data_name))

gg_spearman_power <- ggplot(spearman_power, aes(samples, value, group=sites, col=data, shape=sites, linetype=sites)) +
    ggtitle("Spearman's rank correlation") +
    geom_point() +
    geom_line() +
    scale_y_continuous(name="power", limits=c(0,1)) +
    scale_x_continuous(name="number of samples per site") +
    facet_grid(data~key, labeller=label_parsed) +
    scale_color_manual(values=c("#619CFF", "#00BA38"), guide=FALSE) +
    scale_shape_discrete(name="number of sites") +
    scale_linetype_discrete(name="number of sites") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          panel.grid=element_blank(),
          panel.spacing=grid::unit(0,"lines"))

if (save) ggsave("power.pdf", gg_spearman_power, width=6, height=4)

sumdf <- simdf %>%
    group_by(data, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

sig_sumdf <- simdf %>%
    filter(p.value < level) %>%
    group_by(data, test, sites, samples) %>%
    summarize(
        median=median(effect.size, na.rm=TRUE),
        lwr=quantile(effect.size, 0.025, na.rm=TRUE),
        upr=quantile(effect.size, 0.975, na.rm=TRUE)
    )

(gg_sum <- ggplot(sumdf, aes(samples, median, fill=sites)) +
        geom_ribbon(aes(ymax=upr, ymin=lwr), alpha=0.3) +
        geom_point(aes(col=sites)) +
        geom_line(aes(col=sites)) +
        facet_grid(test~data, scale="free") +
        geom_hline(yintercept=0, lty=2))

gg_sum %+% sig_sumdf

pardf <- clean_list$parlist %>%
    filter(run==3) %>%
    mutate(sim=1:50) %>%
    rename(data=fit)

alldf <- simdf %>% 
    group_by(data,sim,sites,samples, test) %>%
    summarize(
        mean.effect=mean(effect.size, na.rm=TRUE),
        positive.power=mean(p.value<level & effect.size > 0, na.rm=TRUE),
        negative.power=mean(p.value<level & effect.size < 0, na.rm=TRUE)
    ) %>%
    gather(key, value, -data, -sim, -sites, -samples, -test, -mean.effect) %>%
    rename(type=key, power=value) %>%
    merge(pardf) %>%
    filter(samples==100, sites==20)

ggplot(alldf, aes(value, mean.effect)) +
    geom_point(alpha=0.3, aes(col=data)) +
    geom_smooth(method='lm') +
    facet_grid(test~key, scale="free") +
    geom_hline(yintercept=0, lty=2)

pardf_mean <- pardf %>%
    summarize(mean=mean(value))

ggplot(alldf, aes(value, power, col=data)) +
    geom_point(alpha=0.3) +
    geom_smooth(method='lm', aes(fill=data), alpha=0.2) +
    geom_vline(data=pardf_mean, aes(xintercept=mean, col=data), lty=2) + 
    scale_y_continuous(limits=c(0, 1)) +
    facet_grid(test+type~key, scale="free_x")
