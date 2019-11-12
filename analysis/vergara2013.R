library(ggplot2); theme_set(theme_bw())
library(scales) ## for trans_new?
library(dplyr)
library(tidyr)
library(quantreg)

alexandrina <- read.csv("../data/alexandrina.txt", sep="\t")

alexandrina_clean <- alexandrina %>%
    group_by(sites, habitat) %>%
    summarize(male=mean(sex=="male"), 
              infection=mean(Microphallus==1))

alexandrina_clean %>%
    filter(habitat=="shallow") %>%
    ggplot(aes(infection,male)) +
        geom_point() +
        geom_smooth(method="lm")

alexandrina_clean %>%
    group_by() %>%
    filter(habitat=="shallow") %>%
    filter(sites %in% c("campe", "halfway", "wpoint", "swamp")) %>%
    summarize(m=mean(male), i=mean(infection))
