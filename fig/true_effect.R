library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw(base_size = 14,
                                     base_family = "Times"))
library(readxl)
source("../R/powerfun.R")

load("../data/SMC_summary.rda")

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2005)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

logit <- transfun("logit")
arcsin <- transfun("arcsin")

transfun2 <- function(data) {
    list(
        raw=data,
        logit=logit(data),
        arcsin=arcsin(data)
    )
}

test_all <- function(data) {
    data.frame(
        linear=test_lm(data)$effect.size,
        spearman=test_spearman(data)$effect.size,
        quad=test_quad(data)$effect.size,
        quantile_quad=test_quad_rq(data)$effect.size
    )
}

gen <- 1001:1100

simlist <- list(
    dagan=simlist$dagan[[3]],
    vergara=simlist$vergara[[3]],
    mckone=simlist$mckone[[3]]
)

nsim <- length(simlist$dagan)

reslist <- dflist <- vector('list', length(simlist))
names(dflist) <- names(reslist) <- names(simlist)

set.seed(101)
for (sim_name in names(simlist)) {
    siml <- simlist[[sim_name]]
    
    sub_dflist <- sub_reslist <- vector('list', nsim)
    
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
        
        transdf <- transfun2(df)
        
        tres <- lapply(transdf, test_all)
        sub_reslist[[i]] <- tres %>% bind_rows(.id="transformation")
    }
    dflist[[sim_name]] <- sub_dflist
    reslist[[sim_name]] <- sub_reslist %>% bind_rows(.id="sim")
}

resdf <-  reslist %>% 
    bind_rows(.id="data") %>%
    gather(key, value, -data, -sim, -transformation) %>%
    rename(test=key, effect=value) %>%
    as.tbl %>%
    mutate(tgroup=ifelse(test %in% c("linear", "spearman"), "cor", "quad"))

corres <- resdf %>% 
    filter(tgroup=="cor") %>%
    mutate(test=factor(test, labels=c("Pearson correlation", "Spearman's rho")))

## is spearman better?
corres %>%
    mutate(test=factor(test, labels=c("Pearson", "Spearman"))) %>%
    spread(test, effect) %>%
    filter(transformation=="raw") %>%
    group_by(data) %>%
    summarize(prop=mean(Spearman > 0 & Pearson < 0),
              prop2=mean(Spearman < 0 & Pearson > 0))

g1 <- ggplot(corres) +
    geom_boxplot(aes(transformation, effect, col=test)) +
    geom_hline(yintercept = 0, lty=2) +
    scale_y_continuous(name="correlation coefficient") +
    xlab("") +
    facet_grid(data~test, scale="free") +
    theme(legend.position="none",
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          strip.text.y = element_blank())
    
quadres <- resdf %>%
    filter(tgroup=="quad") %>%
    mutate(test=factor(test, labels=c("Least~squares~regression", "Quantile~regression~(tau==0.9)"))) %>%
    mutate(data=factor(data, labels=data_name))

## is logit bad?
quadres %>%
    mutate(test=factor(test, labels=c("LS", "QR"))) %>%
    spread(transformation, effect) %>%
    group_by(data, test) %>%
    summarize(prop.arc=mean(arcsin<0),
              prop.logit=mean(logit<0),
              prop.raw=mean(raw<0))

## is quantile actually better?
quadres %>%
    mutate(test=factor(test, labels=c("LS", "QR"))) %>%
    spread(test, effect) %>%
    filter(transformation=="raw") %>%
    group_by(data) %>%
    summarize(prop=mean(QR<0 & LS >0),
              prop2=mean(QR>0 & LS < 0))

g2 <- ggplot(quadres) +
    geom_boxplot(aes(transformation, effect, col=test)) +
    geom_hline(yintercept = 0, lty=2) +
    xlab("") +
    ylab("standardized quadratic coefficient") +
    facet_grid(data~test, scale="free", labeller = label_parsed) +
    theme(legend.position="none",
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))

geffect <- arrangeGrob(g1, g2, nrow=1, widths=c(0.9, 1))

if (save) ggsave("true_effect.pdf", geffect, width=8, height=6)

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
    scale_alpha_discrete(range=c(0.1, 0), guide=FALSE) +
    scale_color_discrete(guide=FALSE) +
    theme(
        panel.spacing=grid::unit(0,"lines"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        legend.position=c(0.11, 0.82),
        legend.title = element_blank()
    )

## experimental but probably deprecated and not necessary
if (FALSE) {
    pround <- 0.05
    
    comb_sim2 <- comb_sim %>%
        mutate(pinf=round(pinf/pround)*pround, psex=round(psex/pround)*pround) %>%
        group_by(data, pinf, psex) %>%
        tally %>%
        group_by(data) %>%
        mutate(intensity=2*plogis(n/max(n), location=0, scale=0.1)-1)
    
    gsim_tally <- ggplot(comb_sim2, aes(pinf, psex, fill=data)) +
        geom_raster(aes(alpha=intensity)) +
        geom_point(data=aggr, pch=2, size=2.5, col="black") +
        scale_x_continuous(name="Proportion infected", breaks=seq(0,1,0.2), expand=c(0, 0.07)) +
        scale_y_continuous(name="Proportion sexual", breaks=seq(0,1,0.2)) +
        facet_grid(~data, labeller = label_parsed) +
        theme(
            panel.spacing=grid::unit(0,"lines"),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"),
            panel.grid = element_blank(),
            legend.position="none"
        )
}

if (save) ggsave("simulated_data.pdf", gsim, width=7, height=3)

if (save) save("resdf", "comb_sim", file="true_effect.rda")
