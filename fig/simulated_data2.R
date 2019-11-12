library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(readxl)
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
}

source("../R/HPDRegion.R")

save <- FALSE

files <- c("SMC_dagan.rda", "SMC_mckone_assemble.rda", "SMC_vergara.rda")
fits <- c("dagan", "mckone", "vergara")

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

target.gen <- 1099:1100

fit_list <- vector('list', length(fits))

for (i in 1:length(fits)) {
    load(paste0("../data/", files[i]))
    
    par_list <- vector('list', 100)
    
    for (j in 1:10) {
        print(paste(i, j, sep=", "))
        batch_num <- j-1
        
        fn_sim <- paste0("../simdata/post_", fits[i], "_sim_", batch_num, ".rda")
        
        load(fn_sim)
        
        for (k in 1:10) {
            simlist <- sim_reslist[[k]]
            
            df_list <- vector('list', 10)
            
            for (l in 1:10) {
                sim <- simlist[[l]]
                
                S <- sim$S.count[target.gen,]
                SI <- sim$SI.count[target.gen,]
                A <- sim$A.count[target.gen,]
                AI <- sim$AI.count[target.gen,]
                N <- S+A
                
                df_list[[l]] <- data.frame(
                    infected=colMeans((AI+SI)/N),
                    sexual=colMeans(S/N),
                    sim=l
                )
            }
            
            par_list[[10*(j-1)+k]] <- df_list %>%
                bind_rows %>%
                mutate(param=10*(j-1)+k)
        }
        
    }
    
    weight <- data.frame(
        weight=ww[[4]],
        param=1:100
    )
    
    fit_list[[i]] <- par_list %>% 
        bind_rows %>%
        merge(weight)
}

slist <- fit_list %>%
    lapply(bind_rows)

names(slist) <- data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
                               expression(McKone~italic(et~al.)~"(2016)"), 
                               expression(Vergara~italic(et~al.)~"(2014)"))

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

cor(vergara_obs$pinf, vergara_obs$psex, method="spearman")
cor(asin(sqrt(vergara_obs$pinf)), vergara_obs$psex)
cor(dagan$pinf, dagan$psex, method="spearman")
cor(asin(sqrt(dagan$pinf)), asin(sqrt(dagan$psex)))

aggr <- list(mckone=mckone, dagan=dagan, vergara=vergara) %>%
    bind_rows(.id="fit") %>%
    mutate(fit=factor(fit, labels=data_name), sim=NA) 

gg_df <- slist %>%
    bind_rows(.id="fit")

gg_density <- gg_df %>%
    mutate(
        xregion=floor(infected/0.05)*0.05,
        yregion=floor(sexual/0.05)*0.05
    ) %>%
    group_by(fit, xregion, yregion) %>%
    summarize(density=sum(weight)) %>%
    group_by(fit) %>%
    mutate(reldensity=log(density/max(density)))

gsim <- ggplot(gg_density) +
    geom_raster(aes(xregion, yregion, fill=reldensity), hjust = 0, vjust = 0) +
    geom_point(data=aggr,aes(pinf, psex), pch=2, size=2.5, col="black") +
    geom_contour(aes(xregion, yregion, z=reldensity), breaks=-4, col="black", 
                 lty=2) +
    facet_wrap(~fit) +
    scale_x_continuous("Proportion infected", limits=c(0, 1), breaks=seq(0,1,0.2), expand=c(0, 0.07)) +
    scale_y_continuous("Proportion sexual", limits=c(-0, 1), expand=c(0, 0.02)) +
    scale_fill_gradient2(name="log(density)", low="white", mid="#c5c5e8", high="darkred", midpoint=-5) +
    scale_alpha_continuous(guide=FALSE) +
    scale_linetype("probability") +
    facet_wrap(~fit, labeller = label_parsed) +
    theme(
        panel.spacing=grid::unit(0,"lines"),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank()
    )

if (save) ggsave("simulated_data2.pdf", gsim, width=8, height=3)
