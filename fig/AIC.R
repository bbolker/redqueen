library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2); theme_set(theme_bw())
source("../R/powerfun.R")

load("../data/SMC_summary.rda")

save <- FALSE

gen <- 1001:1100

simlist <- list(
    dagan=simlist$dagan[[3]],
    vergara=simlist$vergara[[3]],
    mckone=simlist$mckone[[3]]
)

nsim <- length(simlist$dagan)

reslist <- vector('list', length(simlist))
names(reslist) <- names(simlist)

for (sim_name in names(simlist)) {
    siml <- simlist[[sim_name]]
    
    AICvec <- rep(NA, nsim)
    
    for (i in 1:nsim) {
        sim <- siml[[i]]
        
        df <- with(sim,{
            N.count <- S.count + A.count
            S <- (S.count/N.count)[gen,]
            I <- ((SI.count + AI.count)/N.count)[gen,]
            data.frame(
                infected=colMeans(I),
                sexual=colMeans(S)
            )
        })
        
        l <- lm(sexual~infected, data=df)
        
        l2 <- lm(sexual~I(infected)+I(infected^2), data=df)
        
        AICvec[i] <- AIC(l) - AIC(l2)
    }
    
    reslist[[sim_name]] <-data.frame(AIC=AICvec, sim=1:nsim)
}

resdf <-  reslist %>% 
    bind_rows(.id="data") %>%
    gather(key, value, -data, -sim)

ggplot(resdf) +
    geom_violin(aes(data, value)) +
    geom_hline(yintercept=0, lty=2)

if (save) ggsave("true_effect.pdf", geffect, width=8, height=6)
