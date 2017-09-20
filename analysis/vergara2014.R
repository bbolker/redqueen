library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(MASS) ## for contr.sdif
## load MASS *first*, select() masks dplyr!
library(dplyr)
library(tidyr)
library(car)
library(readxl)

vergara <- read_excel("../data/vergara2014.xlsx")

cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara <- filter(vergara, Ploidy != "UNKNOWN")
vergara_nomale <- vergara %>%
    filter(Gender != "male")

fit <- glm(Microphallus~(Year+Site+Ploidy)^2, 
           family=binomial("logit"), 
           data=vergara_nomale,
           contrasts=list(Year=MASS::contr.sdif,
                          Site=MASS::contr.sdif,
                          Ploidy=MASS::contr.sdif))
Anova(fit, test.statistic = "Wald", type="3")

sum_df <- vergara %>%
    group_by(Year, Site) %>%
    summarize(sex=mean(Ploidy=="sexual"),
              inf=mean(Microphallus)) 

sumfun <- function(x, ...) {
    name <- as.character(substitute(...))
    x %>%
        group_by(...) %>%
        summarize(inf=mean(inf), sex=mean(sex)) %>%
        select(c(inf, sex)) %>%
        summarize_all(funs(mean, sd(.)/mean(.))) %>%
        setNames(c("pinf.mean", "psex.mean", paste0("pinf.", name, "CV"), paste0("psex.", name, "CV"))) %>%
        unlist
}

year_summ <- sumfun(sum_df, Year)
site_summ <- sumfun(sum_df, Site)

vergara_summ <- c(year_summ, site_summ)[c("pinf.YearCV", "psex.YearCV", "pinf.SiteCV", "psex.SiteCV", "pinf.mean", "psex.mean")]

names(vergara_summ) <- c("pinf.timeCV", "psex.timeCV", "pinf.siteCV", "psex.siteCV", "pinf.mean", "psex.mean")

save("vergara_summ", file="vergara_summ.rda")
