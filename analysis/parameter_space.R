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
    res[[i]] <- stochastic_discrete_model(beta=beta_vec[i], bU=10, bI=1.5, seed=103)[c("S.count", "SI.count", "A.count", "AI.count")]
    summary[[i]] <- res[[i]] %>%
        bind_cols() %>%
        filter(row_number() %in% t_sub) %>%
        group_by(id=row_number()) %>%
        summarize(asexual.ratio=ifelse(S.count+A.count==0, 0, A.count/(S.count+A.count)),
                  sexual.infected=ifelse(S.count==0,0, SI.count/S.count),
                  asexual.infected=ifelse(A.count==0, 0, AI.count/A.count)) %>%
        summarize_each(funs(mean(.)), asexual.ratio, sexual.infected, asexual.infected)
}

comb_res <- do.call("rbind", summary)
matplot(beta_vec, comb_res, type="l", log="x", xlab="parasite fecundity", ylab="proportion",
        main="Proportion over 100 generations")
matpoints(beta_vec, comb_res, pch=1)
legend(10, 0.2, 
       legend=c("proportion of asexual", "proportion of infected sexual", "proportion of infected asexual"),
       lty=c(1:3),
       col=c(1:3))
