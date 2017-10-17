library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 14,
                                     base_family = "Times"))
library(gridExtra)
load("../data/SMC_summary.rda")

save <- FALSE

gen <- 1001:1100
vergsim <- 2
sim1 <- simlist$vergara[[3]][[vergsim]]

print(comb_smc$vergara$parlist[[3]][vergsim,])

dflist <- vector('list', 30)

for (i in 1:30) {
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

cclass <- list(
    t1=data.frame(site=c(18, 19)),
    t2=data.frame(site=c(11, 13, 14, 21, 23, 25, 26, 28, 4, 5, 7, 9)),
    t3=data.frame(site=c(12, 16, 24, 27, 3, 30, 6, 8)),
    t4=data.frame(site=c(1, 10, 15, 17, 2, 20, 22, 29))
) %>%
    bind_rows(.id="cycle")

sitedf <- df1 %>%
    filter(site %in% c(1, 11, 16, 18), gen < 1025) %>%
    merge(cclass) %>%
    mutate(cycle=factor(cycle, levels=c("t1", "t2", "t3", "t4"))) %>%
    arrange(cycle)

meandf <- df1 %>%
    group_by(site, key) %>%
    summarize(value=mean(value)) %>%
    spread(key, value)

sampledf <- df1 %>%
    filter(gen==1001) %>%
    spread(key, value)

gcycle <- ggplot(sitedf, aes(gen, value, col=cycle, pch=key, linetype=key)) +
    geom_point() +
    geom_line() +
    scale_y_continuous("proportion", limits=c(0,1), expand=c(0, 0.1)) +
    scale_x_continuous("generation", expand=c(0,0)) +
    scale_color_discrete(guide=FALSE) +
    facet_grid(cycle~.) +
    theme(
        strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position=c(0.42,0.9395),
        legend.title = element_blank(),
        legend.direction = "horizontal"
    )

combdf <- rawdf %>%
    select(-gen) %>%
    bind_rows(meandf, .id="type") %>%
    as.tbl %>%
    mutate(type=factor(type, labels=c("observed effect", "expected effect"))) %>%
    merge(cclass)

gcomb <- ggplot(combdf, aes(infected, sexual)) +
    geom_point(aes(col=cycle, alpha=type, size=type)) +
    geom_smooth(method="lm", formula=y~poly(x, 2), se=FALSE, col=1, lty=2, fullrange=TRUE) +
    scale_x_continuous("proportion infected", limits=c(0, 1), expand=c(0,0)) +
    scale_y_continuous("proportion sexual", limits=c(0, 1)) +
    scale_alpha_discrete(range=c(0.1, 1)) +
    scale_size_discrete(range=c(1, 2.5)) +
    facet_grid(type~.) +
    theme(
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position="none"
    )
    
ggarranged <-arrangeGrob(gcycle, gcomb, nrow=1, widths=c(0.4,0.6))

if (save) ggsave("cycle_example.pdf", ggarranged, width=7, height=4.3)

