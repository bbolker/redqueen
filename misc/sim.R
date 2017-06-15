res <- discrete_model(beta=50, seed=113, migrate="deterministic")

with(res,{
    plot(A.count, type="l", xlim=c(950, 2100))
    lines(S.count, col=2)
})

## proportion of asexual
with(res,{
    N <- A.count+AI.count + S.count+SI.count
    plot((AI.count+A.count)/N, type="l", xlim=c(950, 2100))
})

## prevalence
with(res,{
    plot(AI.count/(A.count+AI.count), type="l", xlim=c(950, 1100))
    lines(SI.count/(S.count+SI.count), col=2)
})

## relative ratio
with(res,{
    u.sexual <- S.count/(S.count+SI.count)
    u.asexual <- A.count/(A.count+AI.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(1000,2100))
    abline(h=1)
})

## host genotype vs parasite genotype
with(res, {
    plot(S[,2,2], type="l", xlim=c(1000, 2000))
    lines(SI[,2,2], col=2)
})

plot(res$SI[,2,2], type="l")
res$S[1100,,]
res$A[1100,,]
res$I[2000,,]

res_spatial <- spatial_discrete_model(beta=50, seed=113, migrate="deterministic", tmax=40000)

## proportion of asexual
with(res_spatial,{
    N <- rowSums(A.count+AI.count + S.count+SI.count)
    plot(rowSums(AI.count+A.count)/N, type="l", xlim=c(950, 2100))
})


## prevalence
with(res_spatial,{
    plot(rowSums(AI.count)/rowSums(A.count+AI.count), type="l", xlim=c(950, 20000))
    lines(rowSums(SI.count)/rowSums(S.count+SI.count), col=2)
})

## relative ratio
with(res_spatial,{
    u.sexual <- rowSums(S.count)/rowSums(S.count+SI.count)
    u.asexual <- rowSums(A.count)/rowSums(A.count+AI.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(39000,40000))
    abline(v=39200, col=2)
    abline(v=39300, col=2)
    abline(v=39400, col=2)
})

