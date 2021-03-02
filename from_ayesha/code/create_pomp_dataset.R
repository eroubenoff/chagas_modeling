# create dataset for pomp
library(lfe)
library(dplyr)
library(lubridate)
library(MASS)
library(Hmisc)
library(pomp)

rm(list = ls())

#### region codes ####
load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata") 

#### cases, births, pop ####
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/births_region_tsir.Rdata") # 1994 - 2016
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/pop_region_tsir.Rdata") # 1992-2016
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_region_tsir.Rdata") #1998 - 2017

ps <- ps %>% mutate(year = as.numeric(year))


tmp <- dat %>% mutate(year = as.numeric(as.character(year)), month = as.character(month)) %>%
  left_join(., all.births, by = c("region", "year", "month")) %>%
  left_join(., ps, by = c("region", "year")) %>%
  filter(as.numeric(as.character(year)) >= 1998 & as.numeric(as.character(year)) <= 2016)


ggplot(data = tmp) + geom_line(aes(x = (year + (as.numeric(month)/12) - (1/12)), y = den/ps, group = region, color = as.factor(region))) + theme_minimal()

#### climate variables ####
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/clim_all.Rdata")
clim <- all_var %>% left_join(., region_master, by = c("code")) %>%
  dplyr::select(region, year, month, tmp, pre) %>%
  group_by(region, year, month) %>%
  dplyr::summarize(tmp = mean(tmp, na.rm = T),
            pre = mean(pre, na.rm = T)) %>%
  filter(year >= 1998 & year <= 2016) %>%
  mutate(month = as.character(month)) %>% ungroup()

tmp2 <- tmp %>% left_join(.,clim, by = c("region", "year", "month"))

#### urbanization ####
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_interp_1996-2017.Rdata")
keycol <- "year"
valuecol <- "perc.urb"
gathercols <- colnames(urb)[2:23]
urb.wide <- gather_(urb, keycol, valuecol, gathercols, na.rm = FALSE) %>%
  mutate(year = as.numeric(as.character(year)))
  
tmp3 <- tmp2 %>% left_join(.,urb.wide, by = c("region", "year"))


#### add splines, time variable and save file for each region ####
for (r in 1:5) {
  data <- tmp3 %>% filter(region == r) %>%
    arrange(year, month) %>%
    mutate(time = year + as.numeric(month)/12 - 1/12) %>%
    dplyr::select(time, den, births, ps, tmp, pre, perc.urb) %>%
    rename(cases = den, pop = ps) %>% arrange(time)
  
  # add splines #
  splinesData <- periodic.bspline.basis(x=data$time,period=1,degree=3,nbasis=6,name="sp%d")
  data <- cbind(data, splinesData)
  
  
  
  save(data, file = paste0("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/forPomp_region_",r,".Rdata"))
}


rm(list = ls())

#### load data ####
r = 3 # central west
load(paste0("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/forPomp_region_",r,".Rdata"))
head(data)
