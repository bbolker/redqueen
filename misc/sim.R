res <- discrete_model(beta=50, tmax=2100, seed=113)

with(res,{
    plot(A.count, type="l", xlim=c(950, 1200))
    lines(S.count, col=2)
})

with(res,{
    plot(AI.count/(A.count+AI.count), type="l", xlim=c(950, 1100))
    lines(SI.count/(S.count+SI.count), col=2)
})

with(res,{
    u.sexual <- S.count/(S.count+SI.count)
    u.asexual <- A.count/(A.count+AI.count)
    plot(u.sexual/u.asexual, type="l", xlim=c(1500,1600))
    abline(h=1)
})

res$S[2000,,]
res$A[1100,,]
res$I[2000,,]
