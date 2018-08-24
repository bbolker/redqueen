library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(gridExtra)

scale_colour_discrete <- function(...,palette="Set1") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Set1") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
}

save <- FALSE

dir <- "../powerdata/"



for (i in 0:9) {
    fn <- paste0("power_vergara_", i, ".rda")
    load(paste0(dir, fn))
    
    
}




