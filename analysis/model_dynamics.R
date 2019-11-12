library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
source("../R/util.R")
source("../R/stochastic_model.R")

target.gen <- 1001:1100

beta_vec <- c(1, 5, 10, 25, 50)

reslist <- vector('list', length(beta_vec))

names(reslist) <- beta_vec

n.site <- 1

CV <- 1

q <- (2*(1:n.site)-1)/(2*n.site)

set.seed(101)
for (beta in beta_vec) {
    sigma <- sqrt(log(CV+1))
    
    meanlog <- log(beta)-1/2*sigma^2
    
    beta2 <- qlnorm(q, meanlog=meanlog, sdlog=sigma)
    
    sim <- stochastic_spatial_discrete_lim_model(
        n.site=n.site,
        n.genotype=1,
        epsilon.site=0.1,
        beta=beta2,
        bI=0.6,
        c_b=1
    )
    ## probably slower but easier this way
    
    sumlist <- vector('list', n.site)
    
    for (i in 1:n.site) {
        sumlist[[i]] <- with(sim,{
            N.count <- S.count+A.count
            I.count <- SI.count+AI.count
            data.frame(
                gen=target.gen,
                psex=(S.count/N.count)[target.gen,i],
                pinf=(I.count/N.count)[target.gen,i]
            )
        })
    }
    reslist[[beta]] <- sumlist %>%
        bind_rows(.id="site")
    
}

resdf <- reslist %>%
    bind_rows(.id="beta") %>%
    mutate(beta=factor(beta, levels=beta_vec)) %>%
    gather(key, value, -beta, -site, -gen) %>%
    as.tbl

ggplot(resdf) +
    geom_line(aes(gen, value, group=interaction(site, key), col=key), lwd=1) +
    scale_x_continuous(limits=c(1001, 1100)) +
    facet_grid(~beta)

