library(tidyr)
library(mgcv)
library(splines)
library(dplyr)
library(ggplot2)
library(rgdal)
library(maptools)
require(descr)
require(RColorBrewer)
require(plotrix)
library(lfe)

rm(list = ls())

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-notificationregion")
filenames <- list.files(pattern = ".Rdata")
yrs.d <- unlist(lapply(filenames, function(x) substr(x,13,16)))

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-notificationregion")
filenames <- list.files(pattern = ".Rdata")
yrs.m <- unlist(lapply(filenames, function(x) substr(x,14,17)))

years <- yrs.d

dat_all <- NA
for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-notificationregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  months <- colnames(data)[4:15]
  months <- substr(months, 1, 3)
  deng <- data
  
  if (years[y] != "2007"){
    load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-notificationregion/MalariaSINAN_",years[y],".Rdata", sep = ""))
    mal <- data
  } else {
    mal <- NA
  }
  
  dat_long <- NA
  
  time <- seq(as.numeric(years[y]), (as.numeric(years[y])+1), 1/12)[1:12]
  
  for (m in 1:12) {
    mth <- paste(months[m], "_d",sep = "")
    tmp <- deng[c("code", mth)]
    
    if (years[y] != "2007"){
      mth <- paste(months[m], "_m",sep = "")
      tmp2 <-  mal[c("code", mth)]
    } else {
      mth <- paste(months[m], "_m",sep = "")
      tmp2 <- as.data.frame(cbind(tmp$code, rep(NA, length = length(tmp$code))))
      names(tmp2) <- c("code", mth)
    }
    joined <- left_join(tmp, tmp2, by = "code")
    joined[is.na(joined)] <- 0
    names(joined)[2:3] <- c("den", "mal")
    
    joined$month = m
    joined$time = time[m]
    
    dat_long <- rbind(dat_long, joined)
    
  }
  dat_long <- dat_long[-1, ]
  dat_long$year <- years[y]
  
  
  dat_all <- rbind(dat_all, dat_long)
  print(y)
  
}

dat_all <- dat_all[-1, ]

# now merge all files
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/clim_all.Rdata")

dat_all$year <- as.numeric(dat_all$year)
tmp1 <- left_join(dat_all, all_var, by = c("code", "year", "month"))

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modis_all.Rdata")
mod.all$NDVI <- as.numeric(mod.all$NDVI)
mod.all$code <- as.numeric(mod.all$code)
mod.all$year <- as.numeric(mod.all$year)

tmp2 <- left_join(tmp1, mod.all, by = c("code", "year", "month"))

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/pop_all.Rdata")
pop.all$year <- as.numeric(pop.all$year)
tmp3 <- left_join(tmp2, pop.all, by = c("code", "year"))
tmp3$ps <- as.numeric(as.character(tmp3$ps))

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")
tmp4 <- left_join(tmp3, region_master, by = "code")


## create lagged variables
tmp5 <- arrange(tmp4, code, year, month)
tmp6 <- tmp5 %>% group_by(code) %>% mutate(lag1_den=lag(den), lag1_mal = lag(mal), lag1_pre = lag(pre), lag1_tmp = lag(tmp))


#### lat and long
#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
municBR@data$order <- 1:nrow(municBR@data)
longlat <- as.data.frame(cbind(municBR@data$code , coordinates(municBR)))
names(longlat) <- c("code", "long", "lat")
longlat$code <- as.numeric(as.character(longlat$code))
tmp7 <- left_join(tmp6, longlat, by = "code")
tmp7$lat <- as.numeric(as.character(tmp7$lat))

# gdp and agric
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/gdp_all.Rdata")
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/agric_all.Rdata")

tmp7 <- left_join(tmp7, gdp.all, by = c("code", "year"))
tmp7 <- left_join(tmp7, agric.all, by = c("code", "year"))


tmp7$month <- as.factor(tmp7$month)
tmp7$year <- as.factor(tmp7$year)
tmp7$stateID <- as.factor(tmp7$stateID)
tmp7$region <- as.factor(tmp7$region)

tmp7$den_inc <- (tmp7$den / tmp7$ps) * 100000
tmp7$mal_inc <- (tmp7$mal / tmp7$ps) * 100000


# create exposure variable
r = tmp7 %>% group_by(month) %>% summarise(r = sum(den, na.rm = TRUE)/sum(ps, na.rm = TRUE))
tmp7 <- left_join(tmp7, r, by = "month")
tmp7$exp.count <- tmp7$ps * tmp7$r

save(tmp7, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic_NOTIFICATION.Rdata")

# plot cases
ts <- tmp7 %>% group_by(region, year, month) %>% summarize(den = sum(den), ps = sum(ps, na.rm = T)) %>%
  mutate(den_inc = den*100000/ps)

ts <- ts[order(ts$region, ts$year, ts$month), ]
head(ts)

ts$time <- rep(seq(2001, 2013, by = 1/12)[-length(seq(2001, 2013, by = 1/12))], 5)

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/timeseries_by_NOTIFICATIONregion.pdf", width=8,height=6)
ggplot(ts, aes(x = time, y = den_inc, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_classic() +
  labs(x = "Month", y = "Monthly dengue incidence per 100,000",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())
dev.off()


ts$time <- rep(seq(2001, 2013, by = 1/12)[-length(seq(2001, 2013, by = 1/12))], 5)

ggplot(ts[ts$time < 2009, ], aes(x = time, y = den/1000, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_classic() +
  labs(x = "Month", y = "Monthly dengue cases (1000)",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())
