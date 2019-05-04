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

data_name <- c(expression(Dagan~italic(et~al.)~"(2013)"), 
               expression(McKone~italic(et~al.)~"(2016)"), 
               expression(Vergara~italic(et~al.)~"(2014)"))

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
    mutate(weight=weight/10,
           plot="violin") 

truesumm <- truedf %>%
    group_by(fit) %>%
    summarize(
        mean=weighted.mean(correlation, weight),
        lwr=wquant(correlation, weight, 0.025),
        upr=wquant(correlation, weight, 0.975)
    )

observed_name <- c(
    expression(McKone~italic(et ~ al.)~"(2016)"),
    expression(Lively~and~Jokela~"(2002)"),
    expression(Kumpulainen~italic(et~al.)~"(2004)"),
    expression(King~italic(et~al.)~"(2011)"),
    expression(Vergara~italic(et~al.)~"(2013)"),
    expression(Gibson~italic(et~al.)~"(2016)")
)

observed_list <- list(
    data.frame(
        value=0.593,
        type="Spearman"
    ),
    data.frame(
        value=0.542,
        type="Pearson"
    ),
    data.frame(
        value=c(0.802, 0.542, 0.811),
        type="Spearman"
    ),
    data.frame(
        value=sqrt(0.37),
        type="Pearson"
    ),
    data.frame(
        value=0.806,
        type="Spearman"
    ),
    data.frame(
        value=c(0.421, 0.617),
        type="Spearman"
    )
)

names(observed_list) <- observed_name

observed <- observed_list %>%
    bind_rows(.id="fit") %>%
    mutate(plot="point")

level <- c(observed_name[-1], data_name)

total <- truedf %>%
    bind_rows(observed) %>%
    mutate(fit=factor(fit, levels=level))

geff <- ggplot(total) +
    geom_violin(aes(fit, correlation, weight=weight, fill=fit), alpha=0.7) +
    geom_point(aes(fit, value, shape=type), size=2) +
    coord_flip() +
    geom_hline(yintercept=0, lty=2) +
    scale_fill_discrete(guide=FALSE) +
    scale_x_discrete("Study", labels=level) +
    scale_y_continuous("Correlation coefficient") +
    scale_shape_manual(values=c(1, 2), label=c("Pearson", "Spearman", "")) +
    theme(
        legend.title = element_blank()
    )

if (save) ggsave("effect_size.pdf", geff, width=6, height=4)
