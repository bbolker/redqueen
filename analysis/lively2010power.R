source("../R/util.R")
source("../R/model.R")

beta_fecundity <- function(beta_vec, model, plot=TRUE, ...) {
    res <- vector("list", length(beta_vec))
    summary <- vector("list", length(beta_vec))
    for(i in 1:length(beta_vec)) {
        cat(i)
        res[[i]] <- model(beta=beta_vec[i], ...)
        summary[[i]] <- apply_gen(res[[i]], c("S.count", "A.count"), mean, gen=c(2001,2100))
    }
    
    if (plot) {
        comb_res <- do.call("rbind", summary)
        matplot(beta_vec, comb_res, type="l", log="xy", xlab="parasite fecundity", ylab="number")
        matpoints(beta_vec, comb_res, pch=1)
    }
    
    return(list(beta_vec=beta_vec, res=res, summary=comb_res))
}

beta_vec <- exp(seq(log(0.1), log(1000), length.out=20))

lively_res <- beta_fecundity(beta_vec, model=lively_model, seed=101)

res <- lively_model(beta=10)

plot(res$S.count, type="l",xlim=c(950,1200))
lines(res$A.count, col=2)

with(res,{
    
})


