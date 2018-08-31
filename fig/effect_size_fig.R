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

source("../R/HPDRegion.R")

save <- FALSE

files <- c("SMC_dagan.rda", "SMC_mckone_assemble.rda", "SMC_vergara.rda")
fits <- c("dagan", "mckone", "vergara")

data_name <- c(expression(Dagan~italic(et~al.)~"(2013) "), 
               expression(McKone~italic(et~al.)~"(2016) "), 
               expression(Vergara~italic(et~al.)~"(2014) "))

target.gen <- 1099:1100

fit_list <- vector('list', length(fits))

for (i in 1:length(fits)) {
    load(paste0("../data/", files[i]))
    
    par_list <- vector('list', 100)
    
    for (j in 1:10) {
        print(paste(i, j, sep=", "))
        batch_num <- j-1
        
        fn_sim <- paste0("../simdata/post_", fits[i], "_sim_", batch_num, ".rda")
        
        load(fn_sim)
        
        for (k in 1:10) {
            simlist <- sim_reslist[[k]]
            
            rho_list <- vector('list', 10)
            
            for (l in 1:10) {
                sim <- simlist[[l]]
                
                S <- sim$S.count[target.gen,]
                SI <- sim$SI.count[target.gen,]
                A <- sim$A.count[target.gen,]
                AI <- sim$AI.count[target.gen,]
                N <- S+A
                
                data <- data.frame(
                    infected=colMeans((AI+SI)/N),
                    sexual=colMeans(S/N)
                )
                
                cc <- cor.test(data$infected, data$sexual, method="spearman")
                
                rho_list[[l]] <- data.frame(
                    correlation=cc$estimate,
                    sim=l
                )
            }
            
            par_list[[10*(j-1)+k]] <- rho_list %>%
                bind_rows %>%
                mutate(param=10*(j-1)+k)
        }
        
    }
    
    weight <- data.frame(
        weight=ww[[4]],
        param=1:100
    )
    
    fit_list[[i]] <- par_list %>% 
        bind_rows %>%
        merge(weight)
}

names(fit_list) <- data_name

truedf <- fit_list %>%
    bind_rows(.id="fit") %>%
    mutate(weight=weight/10) 

truesumm <- truedf %>%
    group_by(fit) %>%
    summarize(
        mean=weighted.mean(correlation, weight),
        lwr=wquant(correlation, weight, 0.025),
        upr=wquant(correlation, weight, 0.975)
    )


ggplot(truedf) +
    geom_violin(aes(fit, correlation, weight=weight, fill=fit), alpha=0.7) +
    coord_flip() +
    theme(
        legend.position = "none"
    )



