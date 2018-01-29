library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
}

save <- FALSE

load("../data/fitted_power.rda")
load("../data/SMC_summary.rda")

powerlist <- reslist

load("../data/SMC_vergara.rda")
load("../data/param_power.rda")

data_name <- c(expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

level <- 0.05

simdf <- powerlist %>%
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

afun <- function(label,size=8, x=25, y=Inf) {
    annotate(geom="text",label=label,x=x,y=y,
             ## http://stackoverflow.com/questions/20083700/how-to-have-annotated-text-style-to-inherit-from-theme-set-options
             family= theme_get()$text[["family"]],
             size=size,
             vjust=1.5,hjust=0)
    ## vjust=0.98,hjust=0.02)
}

g1 <- ggplot(spearman_power %>% filter(key=="positive~correlation"), 
             aes(samples, value, group=sites, col=data, shape=sites, linetype=sites)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(name="power to detect a positive correlation", limits=c(0,1)) +
    scale_x_continuous(name="number of samples per site") +
    scale_color_manual(values=c("#7570B3"), guide=FALSE) +
    scale_shape_discrete(name="number of sites") +
    scale_linetype_discrete(name="number of sites") +
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          panel.grid=element_blank(),
          panel.spacing=grid::unit(0,"lines"),
          legend.position = c(0.75, 0.75))

g1.y <- g1 + afun("a")

g2 <- g1 %+% (spearman_power %>% filter(key=="negative~correlation")) +
    ylab("power to detect a negative correlation") +
    theme(legend.position = "none")

g2.y <- g2 + afun("b")

resdf <- reslist %>% 
    bind_rows %>%
    mutate(power=round(power, digits=2))

g3 <- ggplot(resdf, aes(beta.meanlog, beta.sdlog)) + 
    geom_contour(aes(z=power, colour=..level..), col="black", bins=5) +
    geom_point(data=parlist[[4]], shape=1, col="#7570B3", size=2) +
    scale_x_continuous(expression(beta[meanlog]), limits=c(0,2.5), expand=c(0,0)) +
    scale_y_continuous(expression(beta[sdlog]), limits=c(0, 2), expand=c(0,0)) +
    theme(
        panel.grid = element_blank()
    )

g3.y <- direct.label(g3, 
    list("far.from.others.borders",
        hjust=1.15, vjust=1)
    ) +
    afun("c", x=0.1, y=Inf)

gg_spearman_power <- arrangeGrob(g1.y, g2.y, g3.y, nrow=1)

if (save) ggsave("power.pdf", gg_spearman_power, width=10, height=3.5)
