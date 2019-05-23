library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(emdbook)
source("../R/HPDregion.R")
source("../R/ABC_funs.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

param_list <- vector('list', 3)
names(param_list) <- c("dagan", "mckone", "vergara")

src <- c("../data/")
files <- c("SMC_dagan.rda", "SMC_mckone_assemble.rda", "SMC_vergara.rda")

for (i in 1:3) {
    load(paste0(src, files[i]))
    param_list[[i]] <- cbind(parlist[[4]],
              data.frame(weight=ww[[4]]))
}

data_name <- c(expression(Dagan~italic(et~al.)~"(2013) "), 
               expression(McKone~italic(et~al.)~"(2016) "), 
               expression(Vergara~italic(et~al.)~"(2014) "),
               "Prior")

oldname <- c("beta.mean", "beta.cv", "c_b", "epsilon.site", "V", "n.genotype")
newname <- c("Mean~transmission~rate~(beta[mean])", 
             "CV~transmission~rate~(beta[CV])", 
             "Cost~of~sex~scale~(c[b])", 
             "Mixing~proportion~(epsilon[site])", 
             "Virulence~(V)", 
             "Number~of~asexual~genotypes~(G[asex])")

paramdf <- param_list %>%
    bind_rows(.id="fit") %>%
    gather(key, value, -fit, -weight) %>%
    mutate(value=ifelse(key=="n.genotype", round(value), value)) %>%
    mutate(fit=factor(fit, level=c("dagan", "mckone", "vergara", "Prior"), label=data_name)) %>%
    mutate(key=factor(key, level=oldname, label=newname))

prior <- replicate(10000, as.data.frame(rprior()), simplify=FALSE) %>%
    bind_rows %>%
    gather(key, value) %>%
    mutate(fit="Prior") %>%
    mutate(weight=1/10000) %>%
    mutate(fit=factor(fit, level=c("dagan", "mckone", "vergara", "Prior"), label=data_name)) %>%
    mutate(key=factor(key, level=oldname, label=newname))

combdf <- rbind(paramdf, prior)
 
ggpost <- ggplot(combdf) +
    geom_violin(aes(fit, value, weight=weight, fill=fit), alpha=0.7) +
    coord_flip() +
    facet_wrap(~key, scale="free_x", labeller=label_parsed, nrow=2) +
    scale_fill_discrete(label=c(data_name, "Prior")) +
    theme(
        strip.background = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y=element_blank(),
        axis.title=element_blank(),
        axis.ticks.y=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position="bottom",
        legend.title = element_blank()
    )

if (save) ggsave("posterior.pdf", ggpost, width=8, height=4)

2/wquant(param_list$mckone$c_b, param_list$mckone$weight, c(0.025, 0.975)) 

weighted.mean(2/param_list$mckone$c_b, param_list$mckone$weight)

