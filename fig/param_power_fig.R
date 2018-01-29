library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(directlabels)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

load("../data/SMC_vergara.rda")
load("../data/param_power.rda")

resdf <- reslist %>% 
    bind_rows %>%
    mutate(power=round(power, digits=2))

g <- ggplot(resdf, aes(beta.meanlog, beta.sdlog)) + 
    geom_contour(aes(z=power, colour=..level..), col="black", bins=5) +
    geom_point(data=parlist[[4]], shape=1, col="#7570B3", size=2) +
    scale_x_continuous(expression(beta[meanlog]), limits=c(0,2.5), expand=c(0,0)) +
    scale_y_continuous(expression(beta[sdlog]), limits=c(0, 2), expand=c(0,0)) +
    theme(
        panel.grid = element_blank()
    )

gl <- direct.label(g, 
    list("far.from.others.borders",
        hjust=1, vjust=1.5
    ) 
) 

