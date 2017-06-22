res <- discrete_model(beta=4,bU=30,seed=101)

with(res,{
    plot(A.count, type="l", xlim=c(950, 2100))
    lines(S.count, col=2)
})

## proportion of asexual
with(res,{
    plot((A.count)/N.count, type="l", xlim=c(950, 2100))
})

## prevalence
with(res,{
    plot(AI.count/(A.count), type="l", xlim=c(0, 1600), ylim=c(0,1))
    lines(SI.count/(S.count), col=2)
})

## relative ratio
with(res,{
    u.sexual <- 1-SI.count/(S.count)
    u.asexual <- 1-AI.count/(A.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(1000,2100))
    abline(h=1)
})

## host genotype vs parasite genotype
with(res, {
    plot(AI.count, type="l", xlim=c(1500,1600), ylim=c(0, 20000))
    lines(A.count, col=1, lty=2)
})

plot(res$SI[,2,2], type="l")
res$S[1100,,]
res$A[1100,,]
res$I[2000,,]

res_spatial <- spatial_discrete_model(beta=50, seed=113, migrate="deterministic", tmax=4000)

## proportion of asexual
with(res_spatial,{
    N <- rowSums(A.count+AI.count + S.count+SI.count)
    plot(rowSums(AI.count+A.count)/N, type="l", xlim=c(950, 2100))
})

## prevalence
with(res_spatial,{
    plot(rowSums(AI.count)/rowSums(A.count+AI.count), type="l", xlim=c(950, 4000))
    lines(rowSums(SI.count)/rowSums(S.count+SI.count), col=2)
})

## relative ratio
with(res_spatial,{
    u.sexual <- rowSums(S.count)/rowSums(S.count+SI.count)
    u.asexual <- rowSums(A.count)/rowSums(A.count+AI.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(3900,4000))
})

