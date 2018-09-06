library(MASS)
library(dplyr)
library(tidyr)
library(gridExtra)
library(readxl)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(ellipse)
source("../R/ABC_funs.R")
source("../R/HPDregion.R")

save <- FALSE

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

load_names <- c("dagan", "mckone", "vergara")

dir <- "../simdata/"
slist <- clist <- vector('list', 3)
names(slist) <- names(clist) <- load_names

for (n in load_names) {
    ll <- vector('list', 10)
    
    load(paste0("../data/SMC_", n, ".rda"))
    
    for (i in 0:9) {
        fn <- paste0("post_", n, "_data_",i, ".rda")
        load(paste0(dir, fn))
        
        df <- do.call("rbind", simdata_reslist)
        
        df$weight <- ww[[3]][df$param]
        
        ll[[i + 1]] <- df
    }
    
    slist[[n]] <- do.call("rbind", ll)
    
    slist[[n]]$weight <- slist[[n]]$weight/sum(slist[[n]]$weight)
    
    clist[[n]] <- HPDregion(slist[[n]], w=slist[[n]]$weight, 
                            prob=c(0.95)) %>%
        lapply(function(df) data.frame(
            level=df$level,
            x=df$x,
            y=df$y
        )) %>%
        bind_rows(.id="region")
}

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

names(slist) <- names(clist) <- data_name

mckone <- read.csv("../data/mckone2014.csv", header=FALSE)
vergara <- read_excel("../data/vergara2014.xlsx")
dagan <- read_excel("../data/JEB12245-Dryad.xlsx", sheet=1)

mckone <- (mckone/100) %>%
    setNames(c("infection", "male")) %>%
    mutate(psex=2*round(male, 3),
           pinf=round(infection, 3)) %>%
    dplyr::select(pinf, psex)

vergara <- vergara %>%
    group_by(Year, Site) %>%
    summarize(
        pinf=mean(Microphallus),
        psex=mean(Ploidy=="sexual")
    ) %>%
    group_by() 

vergara_obs <- vergara %>%
    dplyr::select(pinf, psex)

vergara_mean <- vergara %>%
    group_by(Site) %>%
    summarize(
        pinf=mean(pinf),
        psex=mean(psex)
    ) %>%
    group_by() %>%
    dplyr::select(pinf, psex) %>%
    mutate(data=paste(data_name[3]))

dagan <- dagan %>%
    mutate(psex=2*Males) %>%
    rename(pinf=Infected) %>%
    dplyr::select(pinf, psex)

gg_df <- slist %>%
    bind_rows(.id="fit")

gg_contour <- clist %>%
    bind_rows(.id="fit") %>%
    group_by(fit) %>%
    mutate(prob=factor(level, labels=c("95%")))

aggr <- list(mckone=mckone, dagan=dagan) %>%
    bind_rows(.id="fit") %>%
    mutate(fit=factor(fit, labels=data_name[-3]), sim=NA) 

vg_ellipse <- vergara %>%
    split(vergara$Site) %>%
    lapply(function(x){
        with(as.list(x), {
            corr <- cov2cor(cov(data.frame(pinf, psex)))
            ell <- ellipse(corr,
                scale=c(sd(pinf), sd(psex)),
                centre=c(mean(pinf), mean(psex)),
                npoints=100)
            as.data.frame(ell)
        })
    }) %>%
    bind_rows(.id="Site")

vsum <- vergara %>%
    gather(key, value, -Site, -Year) %>%
    group_by(key, Site) %>%
    summarize(
        mean=mean(value),
        lwr=min(value),
        upr=max(value)
    ) %>%
    mutate(fit=as.character(data_name[3]))

vmean <- vsum %>%
    dplyr::select(-upr, -lwr) %>%
    spread(key, mean)

vinfrange <- vsum %>%
    filter(key=="pinf") %>%
    mutate(
        x=vsum$mean[1:4],
        y=vsum$mean[5:8]
    )

vsexrange <- vsum %>%
    filter(key=="psex") %>%
    mutate(
        x=vsum$mean[1:4],
        y=vsum$mean[5:8]
    )

gsim <- ggplot(gg_contour) +
    geom_path(aes(x, y, linetype=prob, col=fit, group=interaction(level, region)), lwd=1) +
    geom_point(data=aggr,aes(pinf, psex), pch=2, size=2.5, col="black") +
    geom_errorbarh(data=vinfrange, aes(x=x, y=y, xmin=lwr, xmax=upr), height=0.02, lty=2) +
    geom_errorbar(data=vsexrange, aes(x=x, y=y, ymin=lwr, ymax=upr), width=0.02, lty=2) +
    geom_point(data=vmean, aes(pinf, psex), pch=1, size=2.5, col="black") +
    scale_x_continuous("proportion infected", limits=c(0, 1), breaks=seq(0,1,0.2), expand=c(0, 0.07)) +
    scale_y_continuous("proportion sexual", limits=c(-0, 1), expand=c(0, 0.02)) +
    scale_colour_discrete(guide=FALSE) +
    scale_linetype("probability") +
    facet_wrap(~fit, labeller = label_parsed) +
    theme(
        panel.spacing=grid::unit(0,"lines"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        legend.position="none"
    )

if (save) ggsave("simulated_data.pdf", gsim, width=7, height=3)
