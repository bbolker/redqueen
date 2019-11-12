library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw())
source("../R/powerfun.R")

load("../data/SMC_summary.rda")

ssimlist <- list(
    dagan=simlist$dagan[[3]],
    vergara=simlist$vergara[[3]],
    mckone=simlist$mckone[[3]]
)

effect_list <- vector('list', length(simlist))
names(effect_list) <- names(simlist)

gen <- c(1001:1100)

for (sim_name in names(simlist)) {
    print(sim_name)
    ss <- ssimlist[[sim_name]]
    
    effect_list2 <- vector('list', 50)
    
    for (k in 1:50) {
        print(k)
        sim <- ss[[k]]
        dflist  <- reslist <- vector('list', length(gen))
        
        for (i in 1:length(gen)) { 
            dflist[[i]] <- tmp.df <- with(sim,{
                N.count <- S.count + A.count
                S <- (S.count/N.count)[gen[i],]
                I <- ((SI.count + AI.count)/N.count)[gen[i],]
                data.frame(
                    infected=I,
                    sexual=S
                )
            })
            
            reslist[[i]] <- data.frame(
                spearman=test_spearman(tmp.df)$effect.size,
                quad=test_quad(tmp.df)$effect.size
            )
        }
        
        names(reslist) <- names(dflist) <- gen
        effect_list2[[k]] <-  reslist %>%
            bind_rows(.id="gen") %>%
            mutate(sim=k)
    }
    effect_list[[sim_name]] <- effect_list2 %>%
        bind_rows()
}

save("effect_list", file="temporal_effect.rda")
