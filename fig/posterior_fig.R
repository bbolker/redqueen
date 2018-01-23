library(emdbook)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(ggstance)
library(gridExtra)
load("../data/SMC_summary.rda")

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

oldname <- c("beta.meanlog", "beta.sdlog", "c_b", "epsilon.site", "V", "n.genotype")
newname <- c("beta[meanlog]", "beta[sdlog]", "c[b]", "epsilon[site]", "V", "G[asex]")

vergara_post <- clean_list$parlist %>% filter(fit=="vergara", run==4)

meanlog.x <- seq(0, 2.8, length=200)
sdlog.x <- seq(0.05, 1.9, length=200)
cb.x <- seq(0.4, 2.2, length=200)
epsilon.x <- seq(0, 0.35, length=200)
V.x <- seq(0.4, 1, length=200)
G.x <- c(0:11)

prior <- rbind(
    data.frame(
        key="beta.meanlog",
        x=meanlog.x,
        y=dcauchy(meanlog.x, location=1, scale=2)
    ),
    data.frame(
        key="beta.sdlog",
        x=sdlog.x,
        y=dlnorm(sdlog.x, meanlog=0, sdlog=2)
    ),
    data.frame(
        key="c_b",
        x=cb.x,
        y=dlnorm(cb.x, meanlog=-0.07, sdlog=0.09)
    ),
    data.frame(
        key="epsilon.site",
        x=epsilon.x,
        y=dbeta(epsilon.x, shape1=1, shape2=9)
    ),
    data.frame(
        key="V",
        x=V.x,
        y=dbeta(V.x, shape1=6, shape2=2)
    ),
    data.frame(
        key="n.genotype",
        x=G.x,
        y=dbetabinom(G.x-1, size=9, prob=3/9, theta=5)
    )
)

prior <- prior %>% 
    mutate(key=factor(key, levels=oldname, labels=newname)) %>%
    group_by(key) %>%
    mutate(y=y/max(y)*1)

vergara_post$key <- factor(vergara_post$key, levels=oldname, labels=newname)

vergara_G <- vergara_post %>%
    filter(key=="G[asex]") %>%
    group_by(fit, run, key, value) %>%
    summarize(count=length(value))

gg_post <- ggplot(prior) +
    geom_line(aes(x, y), lty=2) +
    geom_point(data=prior %>% filter(key=="G[asex]"), aes(x, y)) +
    geom_hline(yintercept = 1.2, alpha=0.15) +
    geom_hline(yintercept = 0, alpha=0.15) +
    geom_boxploth(data=vergara_post %>% filter(key!="G[asex]"), aes(value, 1.2), width=0.2, fill= "#c6c5dd") +
    geom_point(data=vergara_G, aes(value, 1.2, size=count)) +
    scale_size_area(max_size=4) +
    scale_y_discrete(limits=c(0, 1.9), expand=c(0,0.05)) +
    facet_wrap(~key, scale="free_x", labeller=label_parsed, nrow=1) +
    theme(
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, units="cm"),
        legend.position = "none"
    )

ggsave("verg_post.pdf", gg_post, width=7, height=3)
