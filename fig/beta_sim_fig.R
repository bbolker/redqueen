library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(readxl)
load("../data/beta_sim.rda")

mckone <- read.csv("../data/mckone2014.csv", header=FALSE)
vergara <- read_excel("../data/vergara2014.xlsx")
dagan <- read_excel("../data/JEB12245-Dryad.xlsx", sheet=1)

mckone <- (mckone/100) %>%
    setNames(c("infection", "male")) %>%
    mutate(psex=2*round(male, 3),
           pinf=round(infection, 3)) %>%
    select(pinf, psex)

vergara <- vergara %>%
    group_by(Year, Site) %>%
    summarize(
        pinf=mean(Microphallus),
        psex=mean(Ploidy=="sexual")
    ) %>%
    group_by() %>%
    select(pinf, psex)

dagan <- dagan %>%
    mutate(psex=2*Males) %>%
    rename(pinf=Infected) %>%
    select(pinf, psex)

aggr <- list(mckone=mckone, vergara=vergara, dagan=dagan) %>%
    bind_rows(.id="data") %>%
    mutate(asex.genotype="1 asexual genotype")

res <- reslist %>%
    bind_rows(.id="asex.genotype") %>%
    mutate(asex.genotype=factor(asex.genotype, labels=c(
        "1 asexual genotype",
        "2 asexual genotypes",
        "3 asexual genotypes",
        "4 asexual genotypes"
    )))

res3 <- res %>%
    filter(cost %in% c(1, 1.1, 1.5, 2, 2.5, 3))

if (FALSE) {
    
(inf_sex <- ggplot(res3, aes(pinf, psex)) +
    geom_point(aes(col=factor(cost)), alpha=0.5) +
    geom_smooth(aes(col=factor(cost), fill=factor(cost)), alpha=0.2) +
    geom_point(data=aggr, aes(shape=data), size=2.5) +
    scale_x_continuous(name="proportion of infected hosts") +
    scale_y_continuous(name="proportion of sexual hosts") +
    scale_fill_discrete(name="cost of sex") +
    scale_colour_discrete(name="cost of sex") +
    scale_shape_manual(name="data source", 
                       labels=c("Dagan et al., 2005", "McKone et al., 2016", "Vergara et al., 2014"),
                       values=c(15, 16, 17)) +
    facet_wrap(~asex.genotype, nrow=1) +
    theme(
        panel.spacing=grid::unit(0,"lines")
    )) 
}

data_name <- c(expression(Dagan~italic(et~al.)~"(2005)"), 
  expression(McKone~italic(et~al.)~"(2016)"), 
  expression(Vergara~italic(et~al.)~"(2014)"))

(inf_sex2 <- ggplot(res3, aes(pinf, psex)) +
    geom_point(aes(shape=factor(cost)), alpha=0.4) +
    geom_smooth(aes(group=factor(cost), linetype=factor(cost)), alpha=0.15, col="black") +
    geom_point(data=aggr, aes(col=data), size=2) +
    facet_wrap(~asex.genotype, nrow=1) + 
    scale_x_continuous(name="proportion of infected hosts", expand=c(0,0)) +
    scale_y_continuous(name="proportion of sexual hosts") +
    scale_linetype_discrete(name="cost of sex") +
    scale_shape_discrete(name="cost of sex") +
    scale_colour_discrete(name="data source", 
                       labels=data_name) +
    theme(
        panel.spacing=grid::unit(0,"lines")
    ))

ggsave("inf_vs_sex.pdf", inf_sex2, width=9, height=5)

