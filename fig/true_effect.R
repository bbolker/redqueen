library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(readxl)
source("../R/powerfun.R")

load("../data/SMC_summary.rda")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

gen <- 1001:1100

simlist <- list(
    dagan=simlist$dagan[[4]],
    vergara=simlist$vergara[[4]],
    mckone=simlist$mckone[[4]]
)

nsim <- 100

dflist <- vector('list', length(simlist))
names(dflist) <- names(simlist)

set.seed(101)
for (sim_name in names(simlist)) {
    siml <- simlist[[sim_name]]
    
    sub_dflist <- vector('list', nsim)
    
    for (i in 1:nsim) {
        sim <- siml[[i]]
        
        sub_dflist[[i]] <-  df <- with(sim,{
            N.count <- S.count + A.count
            S <- (S.count/N.count)[gen,]
            I <- ((SI.count + AI.count)/N.count)[gen,]
            data.frame(
                infected=colMeans(I),
                sexual=colMeans(S)
            )
        })
    }
    dflist[[sim_name]] <- sub_dflist
}

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
    group_by() 

vergara_obs <- vergara %>%
    select(pinf, psex)

vergara_mean <- vergara %>%
    group_by(Site) %>%
    summarize(
        pinf=mean(pinf),
        psex=mean(psex)
    ) %>%
    group_by() %>%
    select(pinf, psex) %>%
    mutate(data=paste(data_name[3]))

dagan <- dagan %>%
    mutate(psex=2*Males) %>%
    rename(pinf=Infected) %>%
    select(pinf, psex)

aggr <- list(mckone=mckone, vergara=vergara_obs, dagan=dagan) %>%
    bind_rows(.id="data") %>%
    mutate(data=factor(data, labels=data_name), sim=NA) 

comb_sim <- dflist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="data") %>%
    rename(pinf=infected, psex=sexual) %>%
    mutate(data=factor(data, labels=data_name)) %>%
    bind_rows(aggr, .id="type") %>%
    mutate(type= factor(type, labels=c("simulated data", "observed data"))) %>%
    as.tbl

gsim <- ggplot(comb_sim, aes(pinf, psex, col=data)) +
    geom_point(aes(shape=type, alpha=type)) + 
    geom_point(data=aggr, pch=2, size=2.5, col="black") +
    scale_x_continuous(name="proportion infected", breaks=seq(0,1,0.2), expand=c(0, 0.07)) +
    scale_y_continuous(name="proportion sexual", breaks=seq(0,1,0.2)) +
    facet_grid(~data, labeller = label_parsed) +
    scale_shape_manual(values=c(16, 2)) +
    scale_alpha_discrete(range=c(0.2, 0), guide=FALSE) +
    scale_colour_discrete(guide=FALSE) +
    theme(
        panel.spacing=grid::unit(0,"lines"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        legend.position=c(0.11, 0.82),
        legend.title = element_blank()
    )

if (save) ggsave("simulated_data.pdf", gsim, width=7, height=3)

if (save) save("resdf", "comb_sim", file="true_effect.rda")
