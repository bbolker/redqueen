library(dplyr)
source("../R/util.R")
source("../R/model.R")

tmax <- 2100

beta_vec <- exp(seq(log(0.1), log(100), length.out=20))
t_sub <- 2001:2100

res <- vector("list", length(beta_vec))
summary <- vector("list", length(beta_vec))
for(i in 1:length(beta_vec)) {
    cat(i)
    res[[i]] <- discrete_model(beta=beta_vec[i], seed=103)[c("S.count", "SI.count", "A.count", "AI.count")]
    summary[[i]] <- res[[i]] %>%
        bind_cols() %>%
        filter(row_number() %in% t_sub) %>%
        summarize(asexual.ratio=sum(A.count)/sum(A.count+S.count),
                  sexual.infected=sum(SI.count)/sum(S.count),
                  aexual.infected=sum(AI.count)/sum(A.count))
}

comb_res <- do.call("rbind", summary)
matplot(beta_vec, comb_res, type="l", log="x", xlab="parasite fecundity", ylab="proportion",
        main="Proportion over 100 generations")
matpoints(beta_vec, comb_res, pch=1)
legend(10, 0.2, 
       legend=c("proportion of asexual", "proportion of infected sexual", "proportion of infected asexual"),
       lty=c(1:3),
       col=c(1:3))
