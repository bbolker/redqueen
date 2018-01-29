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
    filter(test=="spearman", data=="vergara") %>%
    mutate(data=factor(data, labels=data_name[2]))



g1 <- ggplot(spearman_power %>% filter(key=="positive~correlation"), 
             aes(samples, value, group=sites, col=data, shape=sites, linetype=sites)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(name="power", limits=c(0,1)) +
    scale_x_continuous(name="number of samples per site") +
    scale_color_manual(values=c("#7570B3"), guide=FALSE) +
    scale_shape_discrete(name="number of sites") +
    scale_linetype_discrete(name="number of sites") +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          panel.grid=element_blank(),
          panel.spacing=grid::unit(0,"lines"),
          legend.position = c(0.1, 0.8))

g2 <- g1 %+% (spearman_power %>% filter(key=="negative~correlation")) +
    theme(legend.position = "none")




if (save) ggsave("power.pdf", gg_spearman_power, width=8, height=4.5)
