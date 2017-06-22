source("../R/util.R")
source("../R/model.R")

tmax <- 2100

apply_gen <- function(x, name, FUN, 
                      gen=c(2001,2100)) {
    y <- x[name]
    lapply(y, function(z){
        if(is.vector(z)) {
            FUN(z[gen])
        } else if(is.matrix(z)) {
            FUN(z[gen,])
        } else if(is.array (z)) {
            FUN(z[gen,,])
        }
    })
}

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

beta_vec <- exp(seq(log(0.1), log(100), length.out=20))
discrete_res <- beta_fecundity(beta_vec, model=discrete_model, seed=103)

lively_res <- beta_fecundity(beta_vec, model=lively_model, seed=103)
