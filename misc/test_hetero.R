source("../R/util.R")
source("../R/model.R")

res1 <- discrete_model(seed=101)
res2 <- spatial_discrete_model(n.site=1, epsilon.site =0, seed=101)

## same genotype introduced?
all.equal(
    res1$A[1000,,],
    res2$A[1000,,,1]
)

## same changes in the sexual population?
all.equal(
    res1$S.count,
    res2$S.count[,1]
)

## distribution of infected
all.equal(
    res1$I,
    res2$I[,,1]
)
