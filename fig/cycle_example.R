library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)
load("../data/SMC_summary.rda")

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

gen <- 1001:1100
vergsim <- 3
sim1 <- simlist$vergara[[4]][[vergsim]]

dflist <- vector('list', 40)

for (i in 1:40) {
    dflist[[i]] <- with(sim1,{
        N.count <- S.count + A.count
        S <- (S.count/N.count)[gen,i]
        I <- ((SI.count + AI.count)/N.count)[gen,i]
        data.frame(
            gen=gen,
            infected=I,
            sexual=S
        )
    })
}

rawdf <- dflist %>%
    bind_rows(.id="site")

df1 <- rawdf %>%
    gather(key, value, -gen, -site) %>%
    as.tbl

## look at 1, 10, 11, 12
if (FALSE){
    ggplot(df1, aes(gen, value, col=key)) +
        geom_point() +
        geom_line() +
        facet_wrap(~site)
}

sitedf <- df1 %>%
    filter(site==1) %>%
    filter(gen < 1025) 

meandf <- df1 %>%
    group_by(site, key) %>%
    summarize(value=mean(value)) %>%
    spread(key, value)

gcycle <- ggplot(sitedf, aes(gen, value, pch=key, linetype=key)) +
    geom_point() +
    geom_line() +
    scale_y_continuous("proportion", limits=c(0,1), expand=c(0, 0.1)) +
    scale_x_continuous("generation", expand=c(0,0)) +
    scale_colour_manual(values=cbPalette, guide=FALSE) +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position=c(0.8,0.15),
        legend.title = element_blank(),
        legend.direction = "horizontal"
    )

combdf <- rawdf %>% 
    bind_rows(meandf, .id="type") %>%
    as.tbl %>%
    mutate(type=factor(type, labels=c("accounting temporal variation", "suppressing temporal variation")))

gcomb <- ggplot(combdf, aes(infected, sexual)) +
    geom_point(aes(alpha=type, size=type), col="#7570B3") +
    geom_smooth(method="lm", se=FALSE, col=1, lty=2, fullrange=TRUE) +
    scale_x_continuous("proportion infected", limits=c(0, 1)) +
    scale_y_continuous("proportion sexual", limits=c(0, 1)) +
    scale_colour_manual(values=cbPalette, guide=FALSE) +
    scale_alpha_discrete(range=c(0.025, 1)) +
    scale_size_discrete(range=c(1, 2)) +
    facet_grid(~type) +
    theme(
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position="none"
    )

ggarranged <- arrangeGrob(gcomb, gcycle, nrow=2)

if (save) ggsave("cycle_example.pdf", ggarranged, width=6, height=4.3)

