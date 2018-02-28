simdatapath <- "../data/simdata"

sumlist <- simlist <- parlist <- list()

for (i in 0:29) {
    L <- load(file.path(simdatapath,paste0("SMC_SHARCNET_",i,"_v1.rda")))
    sumlist <- c(sumlist, sumlist2)
    simlist <- c(simlist, simlist2)
    parlist <- c(parlist, parlist2)
}

save("sumlist", file="SMC_sumlist.rda")
save("simlist", file="SMC_simlist.rda")
save("parlist", file="SMC_parlist.rda")
