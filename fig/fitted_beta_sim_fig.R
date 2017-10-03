library(ggplot2); theme_set(theme_bw())
library(dplyr)
library(readxl)
load("../data/fitted_beta_sim.rda")

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

data_name <- c(expression(Dagan~italic(et~al.)~"(2005)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

aggr <- list(mckone=mckone, vergara=vergara, dagan=dagan) %>%
    bind_rows(.id="data") %>%
    mutate(data=factor(data, labels=data_name))

res <- reslist %>%
    bind_rows(.id="data") %>%
    mutate(data=factor(data, labels=data_name)) 

(inf_sex2 <- ggplot(res, aes(pinf, psex)) +
        geom_point(aes(shape=factor(cost)), alpha=0.4) +
        geom_smooth(aes(group=factor(cost), linetype=factor(cost)), alpha=0.15, col="black", span=0.9) +
        geom_point(data=aggr, size=2, col=2) +
        facet_wrap(~data, nrow=1, scale="free_x", labeller=label_parsed) + 
        scale_x_continuous(name="proportion of infected hosts") +
        scale_y_continuous(name="proportion of sexual hosts") +
        scale_linetype_discrete(name="cost of sex") +
        scale_shape_discrete(name="cost of sex") +
        scale_colour_discrete(name="data source", 
                              labels=data_name) +
        theme(
            panel.spacing=grid::unit(0,"lines")
        ))

ggsave("fitted_inf_vs_sex.pdf", inf_sex2, width=9, height=5)

