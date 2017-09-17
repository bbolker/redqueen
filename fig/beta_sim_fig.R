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
    bind_rows(.id="data")

res2 <- res %>%
    mutate(gsex=factor(round(psex, 1)))

ggplot(res2, aes(beta, pinf, col=gsex)) +
    geom_point() +
    geom_smooth(aes(fill=gsex), alpha=0.1) +
    scale_x_log10()

res3 <- res %>%
    filter(cost %in% c(1, 1.1, 1.5, 2, 2.5, 3))

(inf_sex <- ggplot(res3, aes(pinf, psex)) +
    geom_point(aes(col=factor(cost)), alpha=0.5) +
    geom_smooth(aes(col=factor(cost), fill=factor(cost)), alpha=0.2) +
    geom_point(data=aggr, aes(shape=data), size=2.5)) 

ggsave("inf.vs.sex.pdf", inf_sex, width=8, height=6)

