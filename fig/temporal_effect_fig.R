library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

load("../data/SMC_summary.rda")
load("../data/temporal_effect.rda")
load("../data/true_effect.rda")

save <- FALSE

data_name <- c(expression(Dagan~italic(et~al.)~"(2005)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

temporal_effect <- effect_list %>%
    bind_rows(.id="data") %>%
    gather(key, value, -data, -gen, -sim) %>%
    as.tbl %>%
    group_by(data, sim, key) %>%
    summarize(temporal.effect=mean(value, na.rm=TRUE)) %>%
    rename(test=key)
     
true_effect <- resdf %>%
    filter(test %in% c("spearman", "quad"), transformation=="raw")  %>%
    mutate(test=factor(test, labels=c("quad", "spearman"))) %>%
    as.tbl %>%
    select(-c(transformation, tgroup)) %>%
    rename(true.effect=effect)

relative_effect <- merge(true_effect, temporal_effect) %>% 
    as.tbl %>%
    group_by %>%
    mutate(relative.effect=temporal.effect/true.effect)

sensitivity <- clean_list$parlist %>%
    filter(run==3) %>%
    mutate(sim=1:n()) %>%
    group_by() %>%
    select(-run) %>%
    rename(param=key, param.value=value, data=fit) %>%
    as.tbl %>%
    merge(relative_effect)

gsens <- sensitivity %>%
    as.tbl %>%
    select(-relative.effect) %>%
    gather(key, value, -data, -param, -param.value,-test, -sim) %>%
    filter(param %in% c("beta.meanlog", "beta.sdlog")) %>%
    filter(test=="spearman") %>%
    mutate(data=factor(data, labels=data_name),
           param=factor(param, labels=c("beta[meanlog]", "beta[sdlog]")),
           key=factor(key, labels=c("observed~effect", "expected~effect")))

ggsens <- ggplot(gsens, aes(param.value, value)) +
    geom_point(aes(col=data)) +
    geom_smooth(span=0.9, col="black", lty=1) +
    geom_hline(yintercept=0, lty=2) +
    scale_colour_manual(name="data source", 
                        labels=data_name,
                        values=c("#E41A1C", "#4DAF4A", "#377EB8")) +
    scale_x_continuous("parameter value", expand=c(0, 0)) +
    scale_y_continuous("spearman rank correlation") +
    facet_grid(key~param, scale="free", labeller=label_parsed) +
    theme(
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        panel.grid=element_blank(),
        panel.spacing=grid::unit(0,"lines"),
        legend.position="top"
    )

if (save) ggsave("sensitivity.pdf", ggsens, width=8, height=6)

