library(dplyr)
source("../R/util.R")
source("../R/model.R")
source("../R/stochastic_model.R")

beta_vec <- exp(seq(log(0.1), log(100), length.out=20))
t_sub <- 1001:1100

res <- vector("list", length(beta_vec))
summary <- vector("list", length(beta_vec))
for(i in 1:length(beta_vec)) {
    cat(i)
    res[[i]] <- stochastic_spatial_discrete_lim_model(beta=rep(beta_vec[i], 4),seed=103)[c("S.count", "SI.count", "A.count", "AI.count")]
    summary[[i]] <- with(res[[i]],{
        S <- S.count[t_sub,]
        A <- A.count[t_sub,]
        SI <- SI.count[t_sub,]
        AI <- AI.count[t_sub,]
        N <- S+A
        data.frame(
            sexual.ratio=mean(S/N),
            sexual.infected=mean(SI/S),
            asexual.infected=mean(AI/A)
        )
    })
}

comb_res <- do.call("rbind", summary)
matplot(beta_vec, comb_res, type="l", log="x", xlab="parasite fecundity", ylab="proportion",
        main="Proportion over 100 generations")
matpoints(beta_vec, comb_res, pch=1)
legend(0.1, 1, 
       legend=c("proportion of sexual", "proportion of infected sexual", "proportion of infected asexual"),
       lty=c(1:3),
       col=c(1:3))
