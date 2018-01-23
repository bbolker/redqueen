library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

load("../data/SMC_vergara.rda")
load("../data/param_power.rda")

resdf <- reslist %>% 
    bind_rows

ggplot(resdf, aes(beta.meanlog, beta.sdlog)) + 
    geom_raster(aes(fill=power)) +
    geom_point(data=parlist[[4]]) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(limits=c(-0.1, 2.1), expand=c(0,0)) +
    scale_fill_continuous(low="white") +
    theme(
        panel.grid = element_blank()
    )
