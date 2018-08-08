library(emdbook)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)
load("../data/SMC_parlist.rda")

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

oldname <- c("beta.meanlog", "beta.sdlog", "c_b", "epsilon.site", "V", "n.genotype")
newname <- c("beta[meanlog]", "beta[sdlog]", "c[b]", "epsilon[site]", "V", "G[asex]")

vergara_post <- do.call("rbind", parlist) %>%
    as.data.frame %>%
    gather(key, value) %>%
    group_by(key)

vergara_post_sum <- vergara_post %>%
    summarize(
        lwr=quantile(value, 0.025),
        upr=quantile(value, 0.975),
        mean=mean(value)
    ) %>%
    mutate(
        type = "posterior"
    )

vergara_post_sum[5,c("lwr", "upr", "mean")] <- NA

prior <- rbind(
    data.frame(
        key="beta.meanlog",
        lwr=qcauchy(0.025, location=1, scale=2),
        upr=qcauchy(0.975, location=1, scale=2)
    ),
    data.frame(
        key="beta.sdlog",
        lwr=qlnorm(0.025, meanlog=0, sdlog=2),
        upr=qlnorm(0.975, meanlog=0, sdlog=2)
    ),
    data.frame(
        key="c_b",
        lwr=qlnorm(0.025, meanlog=-0.07, sdlog=0.09),
        upr=qlnorm(0.975, meanlog=-0.07, sdlog=0.09)
    ),
    data.frame(
        key="epsilon.site",
        lwr=qbeta(0.025, shape1=1, shape2=9),
        upr=qbeta(0.975, shape1=1, shape2=9)
    ),
    data.frame(
        key="V",
        lwr=qbeta(0.025, shape1=6, shape2=2),
        upr=qbeta(0.975, shape1=6, shape2=2)
    )
)

vergara_post_sum$key <- factor(vergara_post_sum$key, levels=oldname, labels=newname)
prior$key <- factor(prior$key, levels=oldname, labels=newname)
prior$type <- "prior"

pardf <- vergara_post_sum %>%
    bind_rows(prior) %>%
    mutate(type=factor(type, levels=c("prior", "posterior")))

vergara_G <- vergara_post %>%
    filter(key=="n.genotype") %>%
    group_by() %>%
    mutate(key="n.genotype") %>%
    group_by(key, value) %>%
    summarize(count=length(value)) %>%
    mutate(type="posterior")

vergara_G$key <- factor(vergara_G$key, levels=oldname, labels = newname)

prior_G <- data.frame(
    value=1:10,
    count=dbetabinom(0:9,size=9,prob=3/9,theta=5)*100,
    type="prior",
    key="n.genotype") %>%
    mutate(key=factor(key, levels=oldname, labels=newname))

gg_post <- ggplot(pardf) +
    geom_point(aes(mean, type)) +
    geom_point(data=vergara_G, aes(value, type, size=count)) +
    geom_point(data=prior_G, aes(value, type, size=count)) +
    geom_errorbarh(aes(y=type, xmin=lwr, xmax=upr), width=0.2) +
    scale_size_area() +
    scale_x_continuous(expand=c(0.1,0)) +
    facet_wrap(~key, scale="free_x", nrow=2, labeller=label_parsed) +
    theme(
        axis.title = element_blank(),
        strip.background = element_blank(),
        panel.spacing = unit(0, units="cm"),
        legend.position = "none"
    )

if (save) ggsave("verg_post.pdf", gg_post, width=6, height=3)
