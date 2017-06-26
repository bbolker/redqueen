res <- discrete_model(beta=5,seed=113, tmax=2000)

## proportion of asexual
with(res,{
    plot((A.count)/N.count, type="l", xlim=c(1950, 2000))
})

with(res,{
    plot(A.count, type="l", xlim=c(950, 2100))
    lines(S.count, col=2)
})

## prevalence
with(res,{
    plot(AI.count/(A.count), type="l", xlim=c(0, 1500), ylim=c(0,1))
    lines(SI.count/(S.count), col=2)
})

## relative ratio
with(res,{
    u.sexual <- 1-SI.count/(S.count)
    u.asexual <- 1-AI.count/(A.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(1000,2100))
    abline(h=1)
})

with(res, {
    plot(AI.count, type="l", xlim=c(1500,1600), ylim=c(0, 5000))
    lines(A.count, col=1, lty=2)
})

res_spatial <- spatial_discrete_model(beta=c(5:15), n.site=10, seed=113, migrate="deterministic", tmax=2000)

## proportion of asexual
with(res_spatial,{
    plot(rowSums(A.count)/rowSums(N.count), type="l", xlim=c(1950, 2000))
})

## prevalence
with(res_spatial,{
    plot(rowSums(AI.count)/rowSums(A.count), type="l", xlim=c(1900, 2100))
    lines(rowSums(SI.count)/rowSums(S.count), col=2)
})

## relative ratio
with(res_spatial,{
    u.sexual <- 1-rowSums(SI.count)/rowSums(S.count)
    u.asexual <- 1-rowSums(AI.count)/rowSums(A.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(1900,2000), ylim=c(0, 10))
    abline(h=1, lty=2)
})

