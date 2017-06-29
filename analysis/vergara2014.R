library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(MASS) ## for contr.sdif
## load MASS *first*, select() masks dplyr!
library(dplyr)
library(tidyr)
library(car)
library(readxl)

geom_errorbar <- function(...) ggplot2::geom_errorbar(..., width=0.1)
geom_point <- function(...) ggplot2::geom_point(..., size=3)
geom_line <- function(...) ggplot2::geom_line(..., size=1.2)

vergara_theme <- theme(
    axis.line = element_line(colour = "black", size=1.2),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank()
)
vergara <- read_excel("../data/vergara2014.xlsx")

cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara_nomale <- vergara %>%
    filter(Gender != "male", Ploidy != "UNKNOWN")

fit <- glm(Microphallus~(Year+Site+Ploidy)^2, 
           family=binomial("logit"), 
           data=vergara_nomale,
           contrasts=list(Year=MASS::contr.sdif))
Anova(fit, test.statistic = "Wald", type="2")

## Figures

## It seems like Vergara et al are using population standard deviation
## instead of sample standard deviation
## See figure 2 West point 2002
sd2 <- function(x) sd(x) * sqrt((length(x)-1)/length(x))


sites <- c("West Point", "Halfway", "Swamp", "Camp")
vergara_base <- vergara_nomale %>%
    mutate(Site=factor(Site, level=sites), Year=as.numeric(as.character(Year)))

vergara_fig1a <- vergara_base %>% 
    group_by(Year, Ploidy) %>%
    summarize(prev=sum(Microphallus)/length(Microphallus),
              se=sd2(Microphallus)/sqrt(length(Microphallus)))

vergara_fig1b <- vergara_base %>%
    group_by(Year) %>%
    summarize(prop=sum(Ploidy=="asexual")/length(Ploidy),
              se=sd(Microphallus)/sqrt(length(Microphallus)))

fig1a_base <- ggplot(vergara_fig1a, aes(Year, prev, col=Ploidy, group=Ploidy))+
    geom_errorbar(aes(ymin=prev-se, ymax=prev+se), col="black") +
    geom_point() + 
    geom_line() +
    vergara_theme +
    scale_color_grey()

fig1a <- fig1a_base +
    scale_y_continuous(limits=c(0, 0.9), breaks=seq(0, 0.9, by=0.1), expand=c(0,0))

fig1b <- ggplot(vergara_fig1b, aes(Year, prop)) +
    geom_errorbar(aes(ymin=prop-se, ymax=prop+se)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(limits=c(0, 0.6), breaks=seq(0, 0.6, by=0.1), expand=c(0,0)) +
    vergara_theme

grid.arrange(fig1a, fig1b)

## Figure 2
vergara_fig2 <- vergara_base %>%
    group_by(Year, Site, Ploidy) %>%
    summarize(prev=sum(Microphallus)/length(Microphallus), 
              se=sd2(Microphallus)/sqrt(length(Site)),
              Site.count=length(Site))

fig1a_base %+% vergara_fig2 +
    facet_wrap(~Site, scales="free") +
    scale_y_continuous(limits=c(0, 1.05), breaks=seq(0, 1, by=0.1), expand=c(0,0))


## Figure 3
vergara_fig3 <- vergara_fig2 %>%
    filter(Site != "West Point") %>%
    select(-c(se,Site.count)) %>%
    spread(Ploidy, prev) %>%
    summarize(ratio=(1-sexual)/(1-asexual))

ggplot(vergara_fig3, aes(Year, ratio, col=Site, shape=Site)) +
    geom_line() +
    geom_point(size=5,fill="white") +
    geom_hline(yintercept=1, lty=2) +
    ## current order is Halfway, Swamp, Camp
    scale_shape_manual(values=c(22,16,21))+
    scale_y_continuous(limits=c(0, 3.5), breaks=seq(0, 3.5, by=0.5),
                       expand=c(0,0)) +
    vergara_theme +
    scale_color_manual(values=gray(c(0.6,0.4,0.0)))


###########











