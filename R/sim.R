source("util.R")
source("model.R")

n.site <- 50
q <- (2*(1:n.site)-1)/(2*n.site)

beta_vec <- qlnorm(q, meanlog = 1.5, sdlog=0.1)

tmax <- 3000

res_spatial <- spatial_discrete_model(beta=beta_vec, n.site=n.site, seed=113, migrate="deterministic", tmax=tmax)

save("res_spatial", file="50sites.rda")
## move to ../data directory