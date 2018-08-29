library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
}

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2013) "), 
               expression(McKone~italic(et~al.)~"(2016) "), 
               expression(Vergara~italic(et~al.)~"(2014) "))

files <- c("SMC_dagan.rda", "SMC_mckone_assemble.rda", "SMC_vergara.rda")
names <- c("dagan", "mckone", "vergara")

finallist <- vector('list', 3)

for (k in 1:3) {
    load(paste0("../data/", files[k]))
    
    dir <- "../powerdata/"
    
    powerlist <- vector('list', 10)
    
    for (i in 0:9) {
        fn <- paste0("power_", names[k], "_", i, ".rda")
        load(paste0(dir, fn))
        
        powerlist[[i+1]] <- reslist %>%
            bind_rows
    }
    
    weight <- data.frame(
        weight=ww[[4]],
        param=1:100
    )
    
    powerdf <- powerlist %>%
        bind_rows %>%
        merge(weight)
    
    powerdf$fit <- names[k]
    
    finallist[[k]] <- powerdf
}

powersumm <- finallist %>%
    bind_rows %>%
    mutate(site=as.character(site)) %>%
    group_by(sample, site, fit) %>%
    summarize(
        power=weighted.mean(p.value<0.05 & effect.size > 0, w=weight, na.rm=TRUE)
    ) %>%
    mutate(fit=factor(fit, label=data_name))

gpower <- ggplot(powersumm) +
    geom_point(aes(sample, power, shape=site, col=fit)) +
    geom_line(aes(sample, power, group=site, lty=site, col=fit)) +
    scale_y_continuous("power", limit=c(0, 1)) +
    scale_x_continuous("number of samples per site") +
    scale_colour_discrete(guide=FALSE) +
    scale_linetype_discrete("number of sites") +
    scale_shape_discrete("number of sites") +
    facet_wrap(~fit, labeller=label_parsed) +
    theme(
        strip.background = element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position=c(0.9, 0.73)
    )

ggsave("power.pdf", gpower, width=7, height=3.5)
