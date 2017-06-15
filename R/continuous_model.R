initialize <- function(p=0.5, q=0.5, 
                       size=1000) {
    locus1 <- c(p, 1-p)
    locus2 <- c(q, 1-q)
    
    gamete <- unlist(lapply(locus1, function(x) x*locus2))
    genotype <- outer(gamete, gamete, "*")
    genotype <- size * scaled_matrix(genotype)
}

yini <- list(
    XA_S=introduce(),
    XA_I=matrix(0, 4, 4),
    XS_S=initialize(size=1000),
    XS_I=initialize(size=10),
    Y_S=500,
    Y_I=initialize(size=10)
)

pars <- c(
    c_b = 1,
    h_b = 1/5000,
    r = 0.5,
    beta_X = 1,
    mu_X = 1/20,
    alpha=1,
    rho_X=0.05,
    a = 1/20,
    beta_Y = 5,
    mu_Y = 1/20,
    rho_Y=0.05,
    f = 1
)

continuous_model <- function(t, yini, parameters) {
    with(as.list(c(yini, parameters)), {
        XA_S <- matrix(yini[1:16], 4, 4)
        XA_I <- matrix(yini[17:32], 4, 4)
        XS_S <- matrix(yini[33:48], 4, 4)
        XS_I <- matrix(yini[49:64], 4, 4)
        Y_I <- matrix(yini[66:81], 4, 4)
        
        N_X <- scaled_sum(XA_S+XA_I+XS_S+XS_I)
        N_Y <- Y_S + scaled_sum(Y_I)
        
        b_A <- r * (XA_S + f * XA_I) * (1- h_b * N_X)
        b_S <- c_b * r * outcross(XS_S + f * XS_I, rho_X)/2  * (1- h_b * N_X)
        
        Z <- outcross(Y_I, rho_Y)
        
        dXA_S <- b_A - beta_X * XA_S * Z/N_X - mu_X * XA_S - a * XA_S * N_Y/N_X
        dXA_I <- beta_X * XA_S * Z/N_X - (mu_X + alpha) * XA_I - a * XA_I * N_Y/N_X
        
        dXS_S <- b_S - beta_X * XS_S * Z/N_X - mu_X * XS_S - a * XS_S * N_Y/N_X
        dXS_I <- beta_X * XS_S * Z/N_X - (mu_X + alpha) * XS_I - a * XS_I * N_Y/N_X
        
        dY_S <- mu_Y * N_Y - beta_Y * Y_S * scaled_sum(XA_I + XS_I)/N_Y - mu_Y * Y_S
        dY_I <- beta_Y * Y_S * (XA_I + XS_I)/N_Y - mu_Y * Y_I
        
        list(c(dXA_S, dXA_I, dXS_S, dXS_I, dY_S, dY_I))
        
    })
}

## for experimental purposes
## move to a different directory
library(deSolve)
r <- ode(unlist(yini), seq(0, 1000), continuous_model, pars)
r <- as.data.frame(r)

plot(r$XA_S11)
plot(r$XS_S12)
plot(r$Y_S)

