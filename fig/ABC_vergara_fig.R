library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2); theme_set(theme_bw())
library(grid)
library(gridExtra)
load("../data/ABC_vergara.rda")
load("../data/vergara_fitted_power.rda")

vergara <- read_excel("../data/vergara2014.xlsx")
cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara <- filter(vergara, Ploidy != "UNKNOWN")
vergara_nomale <- vergara %>%
    filter(Gender != "male")

vergara_CV <- structure(
    c(0.118560462285074, 0.129710635277051, 0.291586199244854, 0.267463972926705),
    .Names = c("pinf.timeCV", "psex.timeCV", "pinf.siteCV", "psex.siteCV"
    ))

vergara_mean <- structure(
    c(0.433838731993719, 0.703887586155451), 
    .Names = c("pinf.mean", 
               "psex.mean"))

vergara_summ <- c(vergara_CV, vergara_mean)

new.tol <- 0.36 ## sqrt(6*0.15^2)

## sum(sapply(res, inherits, "try-error"))

param.keep2 <- sapply(res, function(x) ifelse(is.null(x) || inherits(x, "try-error"), FALSE, sqrt(sum((vergara_summ - x)^2)) < new.tol))

sum(param.keep2)

param.table2 <- cbind(param.keep2, param.table)

prior <- param.table2 %>%
    as_data_frame %>%
    gather(key, value, -param.keep2)

posterior <- prior %>%
    filter(param.keep2==TRUE) %>%
    select(-param.keep2)

ggplot(prior, aes(x=value)) +
    geom_density() +
    geom_density(data=posterior, col=2) +
    facet_wrap(~key, scales="free")

vergara_summ2 <- as.data.frame(vergara_summ) %>%
    rename(value=vergara_summ) %>%
    mutate(key=rownames(.))

res[param.keep2] %>%
    lapply(. %>% t %>% as.data.frame) %>%
    bind_rows(.id="sim") %>%
    gather(key, value, -sim) %>%
    ggplot() +
    geom_density(aes(value, col=key)) +
    geom_vline(data=vergara_summ2, aes(xintercept=value), lty=2) +
    facet_wrap(~key, scale="free_y")
    
posterior %>%
    group_by(key) %>%
    summarise(mean=mean(value))

nsamples <- seq(10, 100, by=10)
nsites <- seq(10, 50, by=10)

power_site <- list(power_list, effect_list) %>%
    lapply(. %>% bind_rows(., .id="nsites")) %>%
    Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by=c("nsites", "nsamples")), .)

power_site$nsites <- factor(power_site$nsites, labels=nsites)

gg_power <- ggplot(power_site, aes(nsamples, power, col=nsites)) +
    geom_point() +
    geom_line() +
    xlab("") +
    theme(legend.position = "none")

gg_effect <- ggplot(power_site, aes(nsamples, median, col=nsites, fill=nsites)) +
    geom_point() +
    geom_line() +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1) +
    xlab("") +
    ylab("effect size")

grid.arrange(gg_power, gg_effect, nrow=1,
    bottom=textGrob("number of samples per site")
)

# plot(gen_power, type="l")

# matplot(gen_effect[,-1], type="l")

summary_df <- vergara_fitted_summary %>%
    lapply(. %>% cbind(data.frame(gen=1001:1100))) %>%
    bind_rows(.id="run") %>%
    as.tbl %>%
    gather(key,value, -run, -gen) %>%
    spread(run, value) %>%
    apply(1, function(x){
        xx <- as.numeric(x[-c(1:2)])
        data.frame(
            gen=x[1],
            type=x[2],
            median=median(xx),
            lwr=quantile(xx, 0.025),
            upr=quantile(xx, 0.975)
        )
    }) %>%
    bind_rows %>%
    mutate(gen=as.numeric(gen))

ggplot(summary_df,aes(gen, median, col=type, fill=type)) +
    geom_line() +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1)


