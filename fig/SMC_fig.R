library(emdbook)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw(base_size = 12,
                                     base_family = "Times"))
library(grid)
library(gridExtra)

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dakr2") scale_fill_brewer(...,palette=palette)

if (.Platform$OS.type=="windows") {
    windowsFonts(Times=windowsFont("Times"))
} 

save <- FALSE

src <- "../data/"
SMC_files <- c("SMC_vergara.rda", "SMC_dagan.rda", "SMC_mckone_assemble.rda")
summ_files <- c("vergara_summ.rda", "dagan_summ.rda", "mckone_summ.rda")

get_all <- function(x) {
    ## cat(x,"\n")
    L <- load(x)
    r <- lapply(L,function(x) get(x,envir=parent.frame(2)))
    names(r) <- L
    return(r)
}

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

data_name_simple <- c("D", "M", "V")

comb_smc <- comb_summ <- vector('list', 2)

comb_smc <- lapply(paste0(src, SMC_files), get_all)
comb_summ <- lapply(paste0(src, summ_files), get_all)

names(comb_smc) <- names(comb_summ) <- c("vergara", "dagan", "mckone")

clean_fun <- function(x, target) {
    x %>%
        lapply("[[", target) %>%
        lapply(bind_rows, .id="run") %>%
        bind_rows(.id="fit") %>%
        as.tbl %>%
        gather(key, value, -run, -fit) %>%
        group_by(run, fit, key)
}  

clean_list <- c("parlist", "sumlist") %>%
    sapply(clean_fun, x=comb_smc, simplify=FALSE) %>%
    lapply(mutate, param=1:100)

simlist <- comb_smc %>%
    lapply("[[", "simlist")

SMC_weight <- comb_smc %>%
    lapply("[[", "ww") %>%
    lapply(function(x) lapply(x, function(y) {data.frame(param=1:100, weight=y) })) %>%
    lapply(bind_rows, .id="run") %>%
    bind_rows(.id="fit") %>%
    as.tbl

weighted_list <- clean_list %>%
    lapply(merge, SMC_weight)

SMC_summary <- weighted_list %>%
    lapply(group_by, run, fit, key) %>%
    lapply(summarize, value=weighted.mean(value, w=weight, na.rm=TRUE)) %>%
    lapply(spread, key, value)

if(save) save("comb_summ", "weighted_list", "SMC_weight", "SMC_summary", file="SMC_summary.rda")
if(save) save("comb_summ", file="comb_summ.rda")
