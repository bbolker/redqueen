library(readxl)
library(dplyr)
vergara <- read_excel("../data/vergara2014.xlsx")

cols <- c("Year", "Site", "Gender", "Ploidy")
vergara[cols] <- lapply(vergara[cols], factor)
vergara_clean <- vergara %>%
    filter(Ploidy != "UNKNOWN")

vergara_p <- vergara_clean %>%
    group_by(Site, Year) %>%
    ## mean prevalence of infection/sex for each site every year
    summarize(pinf=mean(Microphallus), psex=sum(Ploidy=="sexual")/length(Ploidy))

vergara_year_CV <- vergara_p  %>% 
    group_by(Year) %>%
    ## CV across site
    summarize_each(funs(mean(.)), pinf, psex) %>%
    ## mean of CV across year
    select(-Year) %>%
    summarize_all(funs(sd(.)/mean(.))) %>%
    unlist

vergara_site_CV <- vergara_p  %>% 
    group_by(Site) %>%
    ## CV across year
    summarize_each(funs(mean(.)), pinf, psex) %>%
    ## mean of CV across site
    select(-Site) %>%
    summarize_all(funs(sd(.)/mean(.))) %>%
    unlist

vergara_CV <- c(vergara_year_CV, vergara_site_CV)

vergara_mean <- vergara_p %>%
    group_by %>%
    select(-Site, -Year) %>%
    summarize_each(funs(mean(.))) %>%
    unlist
    
    
