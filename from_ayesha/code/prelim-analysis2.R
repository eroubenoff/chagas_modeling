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


# setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
# filenames <- list.files(pattern = ".Rdata")
# yrs.d <- unlist(lapply(filenames, function(x) substr(x,13,16)))
# 
# setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion")
# filenames <- list.files(pattern = ".Rdata")
# yrs.m <- unlist(lapply(filenames, function(x) substr(x,14,17)))
# 
# years <- intersect(yrs.d, yrs.m)
# 
# dat_all <- NA
# for (y in 1:length(years)){
#   load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
#   months <- colnames(data)[4:15]
#   months <- substr(months, 1, 3)
#   deng <- data
#   
#   load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_",years[y],".Rdata", sep = ""))
#   mal <- data
#   dat_long <- NA
#   
#   time <- seq(as.numeric(years[y]), (as.numeric(years[y])+1), 1/12)[1:12]
#   
#   for (m in 1:12) {
#     mth <- paste(months[m], "_d",sep = "")
#     tmp <- deng[c("code", mth)]
#     
#     mth <- paste(months[m], "_m",sep = "")
#     tmp2 <-  mal[c("code", mth)]
#     
#     joined <- left_join(tmp, tmp2, by = "code")
#     joined[is.na(joined)] <- 0
#     names(joined)[2:3] <- c("den", "mal")
#     
#     joined$both <- ifelse(joined[,2] >= 1 & joined[,3] >= 1, 1, 0)
#     joined$month = m
#     joined$time = time[m]
#     
#     dat_long <- rbind(dat_long, joined)
#     
#   }
#   dat_long <- dat_long[-1, ]
#   dat_long$year <- years[y]
# 
#   
#   dat_all <- rbind(dat_all, dat_long)
#   print(y)
#   
# }
# 
# dat_all <- dat_all[-1, ]
# 
# 
# save(dat_all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_malaria_cooccurence.Rdata")



setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
yrs.d <- unlist(lapply(filenames, function(x) substr(x,13,16)))

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
yrs.m <- unlist(lapply(filenames, function(x) substr(x,14,17)))

years <- yrs.d

dat_all <- NA
for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  months <- colnames(data)[4:15]
  months <- substr(months, 1, 3)
  deng <- data
  
  if (years[y] != "2007"){
    load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_",years[y],".Rdata", sep = ""))
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


save(dat_all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_malaria_cooccurence.Rdata")

################ climate data
var = c("pre", "tmn", "tmp", "tmx")
rm(long, tmp)
for (v in 1:length(var)) {
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/climdata/cleaned_",var[v],".Rdata", sep = "")
  load(filepath)
  
  # reshape to long format
  keycol <- "date"
  valuecol <- var[v]
  gathercols <- colnames(data.all)[which(substr(colnames(data.all),1,1) == "Y")]
  long <- gather_(data.all, keycol, valuecol, gathercols)
  long$year <- as.numeric(substr(long$date,2,5))
  long$month <- as.numeric(substr(long$date,7,8))
  long$code <- as.numeric(long$code)
  if (v == 1) tmp <- long
  if(v>= 2) {
    all_var <- left_join(long, tmp, by = c("code_long", "code", "date", "year", "month"))
    tmp <- all_var
  }
}

save(all_var, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/clim_all.Rdata")


################ modis data

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/")
filenames <- list.files(pattern = "cleaned_NDVI_")
years <- unique(unlist(lapply(filenames, function(x) substr(x,14,17))))

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_2001_6.Rdata")
tmp <- data[0,]
rm(data)


for (m in 6:12){
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_",2001,"_",m,".Rdata", sep = "")
  load(filepath)
  data$month = m
  data$year = 2001
  tmp <- rbind(tmp,data)
  rm(data)
  
}


for (y in 2:length(years)){
  for (m in 1:12){
    filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_",years[y],"_",m,".Rdata", sep = "")
    load(filepath)
    data$month = m
    data$year = years[y]
    tmp <- rbind(tmp,data)
    rm(data)
    
  }
}

mod.all <- tmp[-1,]

save(mod.all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modis_all.Rdata")

################ population data
rm(list = ls())
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")
# reshape to long format
keycol <- "year"
valuecol <- "ps"
gathercols <- colnames(pop)[which(substr(colnames(pop),1,1) == "2")]
pop.all <- gather_(pop, keycol, valuecol, gathercols)
pop.all$code <- as.numeric(pop.all$code)
pop.all$ps <- as.numeric(as.character(pop.all$ps))


# calculate density by year
# area by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")
pop.all <- left_join(pop.all, area, by = c("code"))
pop.all$area_km2 <- as.numeric(as.character(pop.all$area_km2))
pop.all$dens <- pop.all$ps/pop.all$area_km2

# urbanization
# load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")
# pop.tmp <- pop[,c("code", "2000", "2010")]
# names(pop.tmp) <- c("code", "pop2000", "pop2010")
# pop.tmp$pop2000 <- as.numeric(as.character((pop.tmp$pop2000)))
# pop.tmp$pop2010 <- as.numeric(as.character((pop.tmp$pop2010)))
# 
# # urban pop by munic 2000
# load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned_1991-2000.Rdata")
# urban2 <- urban[,c("code", "urban2000")]
# tmp2 <- left_join(pop.tmp, urban2, by = c("code"))
# tmp2$urbanPop <- as.numeric(as.character(tmp2$urban2000))
# tmp2$urban2000 <- (tmp2$urbanPop/tmp2$pop2000) * 100
# tmp2$urban2000[tmp2$urban2000 > 100] <- 100
# tmp3 <- tmp2[,c("code", "urban2000")]
# 
# # urbanizaton % by munic 2010
# load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned.Rdata")
# tmp4 <- left_join(tmp3, urban, by = c("code"))
# tmp4$urban2010 <- as.numeric(as.character(tmp4$Urban))
# tmp4 <- tmp4[,c("code", "urban2000", "urban2010")]
# #plot(tmp4$urban2000, tmp4$urban2010)

load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_interp_2000-2010.Rdata")
keycol <- "year"
valuecol <- "urban.perc"
gathercols <- colnames(urb)[which(substr(colnames(urb),1,1) == "2")]
urb.all <- gather_(urb, keycol, valuecol, gathercols)
urb.all$code <- as.numeric(urb.all$code)
urb.all$urban.perc <- as.numeric(as.character(urb.all$urban.perc))


pop.all <- left_join(pop.all, urb.all, by = c("code","code_long", "year"))

save(pop.all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/pop_all.Rdata")


################ gdp per capita and agricultural value
rm(list = ls())
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/gdp_cleaned.Rdata")
# reshape to long format
keycol <- "year"
valuecol <- "gdp"
gathercols <- colnames(gdp)[which(substr(colnames(gdp),1,1) == "X")]
gdp.all <- gather_(gdp, keycol, valuecol, gathercols)
gdp.all$code <- as.numeric(gdp.all$code)
gdp.all$gdp <- as.numeric(as.character(gdp.all$gdp))
gdp.all$year <- as.numeric(substr(gdp.all$year,2,5))
save(gdp.all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/gdp_all.Rdata")


rm(list = ls())
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/agric_cleaned.Rdata")
# reshape to long format
keycol <- "year"
valuecol <- "agric"
gathercols <- colnames(agric)[which(substr(colnames(agric),1,1) == "X")]
agric.all <- gather_(agric, keycol, valuecol, gathercols)
agric.all$code <- as.numeric(agric.all$code)
agric.all$agric <- as.numeric(as.character(agric.all$agric))
agric.all$year <- as.numeric(substr(agric.all$year,2,5))
save(agric.all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/agric_all.Rdata")


################ Merge all the data
rm(list = ls())
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_malaria_cooccurence.Rdata")
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

save(tmp7, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic.Rdata")



############ fit model

reg <- felm(den ~ offset(ps) + NDVI + tmp + pre + dens + log(gdp)  | month:stateID   | 0 | region, data=tmp7)
summary(reg)

reg <- felm(mal ~ offset(ps) + NDVI + tmp + pre + dens + log(gdp) | month:stateID   | 0 | region, data=tmp7)
summary(reg)


reg <- felm(den ~ offset(exp.count) + NDVI + tmp + pre + dens + log(gdp) + urban2000  | month:stateID   | 0 | region, data=tmp7)
summary(reg)


#check spatial autocorrelation - appears to be a spatial structure
tmp8 <- tmp7[,c("den","ps","NDVI","lag1_tmp", "lag1_pre", "month", "stateID", "lat")]
tmp8 <- tmp8[complete.cases(tmp8[,1:7]),]
plot(tmp8$lat, reg$residuals)


fit <- glm(den ~ offset(ps) + NDVI + tmp + pre + month + stateID, data = tmp4, family = "binomial")
summary(fit)

fit <- glm(den ~ NDVI + tmp + pre + dens + log(gdp) + month + stateID, data = tmp7,
           family = poisson(link = "log"),
           offset = log(exp.count))
summary(fit)

##### make plots

### temperature
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_temp.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$tmp, (tmp7$den_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "Temp (C)")
plot(tmp7$tmp, (tmp7$mal_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "Temp (C)")
legend("topleft",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)

dev.off()

### precipitation
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_prec.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$pre, (tmp7$den_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region],pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "Precipitation (mm)")
plot(tmp7$pre, (tmp7$mal_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region],pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "Precipitation (mm)")
legend("topright",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)
dev.off()

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_prec_lag1.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$lag1_pre, (tmp7$den_inc), pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "Lagged Precipitation (mm)")
plot(tmp7$lag1_pre, (tmp7$mal_inc), pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "Lagged Precipitation (mm)")
dev.off()

### vegetation index
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_NDVI.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$NDVI, (tmp7$den_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region],pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "Vegetation index")
plot(tmp7$NDVI, (tmp7$mal_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "Vegetation index")
legend("topleft",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)
dev.off()
#0 = no greenness, 1 = maximum green

### density
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_dens.png", 640, 480)
par(family = "serif")
par(mfrow = c(2,2))
plot(tmp7$dens, (tmp7$den_inc), pch = 20, col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], xlim = c(0,100),main = "Dengue", ylab = "cases per 100,000", xlab = "Density (people/km2)")
plot(tmp7$dens, (tmp7$mal_inc), pch = 20, col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], xlim = c(0,100),main = "Malaria", ylab = "cases per 100,000", xlab = "Density (people/km2)")

plot(tmp7$dens, (tmp7$den_inc), pch = 20, col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], main = "Dengue", ylab = "cases per 100,000", xlab = "Density (people/km2)")
plot(tmp7$dens, (tmp7$mal_inc), pch = 20, col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], main = "Malaria", ylab = "cases per 100,000", xlab = "Density (people/km2)")
legend("topright",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)
dev.off()

### gdp per capita
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_gdp.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$gdp, (tmp7$den_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "GDP per capita")
plot(tmp7$gdp, (tmp7$mal_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "GDP per capita")
legend("topleft",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)

dev.off()

# agricultural value
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_agric.png", 640, 480)
par(family = "serif")
par(mfrow = c(1,2))
plot(tmp7$agric, (tmp7$den_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Dengue", ylab = "cases per 100,000", xlab = "Value from Agricultural")
plot(tmp7$agric, (tmp7$mal_inc), col = c("red", "blue", "yellow", "green", "orange")[tmp7$region], pch = 20, main = "Malaria", ylab = "cases per 100,000", xlab = "Value from Agricultural")
legend("topleft",  c("North", "Northeast", "Central West", "South", "Southeast"), col = c("red", "blue", "yellow", "green", "orange"), pch = 20)

dev.off()

#### separate out regions
### temperature
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_temp_reg.png", width=6,height=4,units="in",res=1200)
par(family = "serif"); par(mar = c(4,3,3,2)); par(mgp = c(2,1,0))
par(mfcol = c(2,5))
for (r in 1:5) {
  plot(tmp7$tmp[tmp7$region == paste(r)], tmp7$den_inc[tmp7$region == paste(r)], xlim = c(5,35), pch = 20, main = paste("Dengue, Region ",r), ylab = "cases per 100,000", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7, xlab = "")
  plot(tmp7$tmp[tmp7$region == paste(r)], tmp7$mal_inc[tmp7$region == paste(r)], xlim = c(5,35), pch = 20, main = paste("Malaria, Region ",r), ylab = "cases per 100,000", xlab = "Temp (C)", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7)
  
}
dev.off()

### gdp per capita
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_gdp_reg.png", width=6,height=4,units="in",res=1200)
par(family = "serif"); par(mar = c(4,3,3,2));  par(mgp = c(2,1,0))
par(mfcol = c(2,5))
for (r in 1:5) {
  plot(tmp7$gdp[tmp7$region == paste(r)], tmp7$den_inc[tmp7$region == paste(r)],  pch = 20, main = paste("Dengue, Region ",r), ylab = "cases per 100,000", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7, xlab = "")
  plot(tmp7$gdp[tmp7$region == paste(r)], tmp7$mal_inc[tmp7$region == paste(r)], pch = 20, main = paste("Malaria, Region ",r), ylab = "cases per 100,000", xlab = "GDP", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7)
  
}
dev.off()

### vegetation index
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_NDVI_reg1.png", width=6,height=4,units="in",res=1200)
par(family = "serif"); par(mar = c(4,3,3,2));  par(mgp = c(2,1,0))
par(mfcol = c(2,5))
for (r in 1:5) {
  plot(tmp7$NDVI[tmp7$region == paste(r)], tmp7$den_inc[tmp7$region == paste(r)], xlim = c(0,1), pch = 20, main = paste("Dengue, Region ",r), ylab = "cases per 100,000", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7, xlab = "")
  plot(tmp7$NDVI[tmp7$region == paste(r)], tmp7$mal_inc[tmp7$region == paste(r)], xlim = c(0,1), pch = 20, main = paste("Malaria, Region ",r), ylab = "cases per 100,000", xlab = "NDVI", cex.lab = 0.7, cex.main = 0.7, cex.axis = 0.7)
  
}
dev.off()

library(ggplot2)
library(cowplot)
p1 <- ggplot(tmp7, aes(x=NDVI, y=den_inc)) + geom_point(size = 0.5, alpha = 1/5) + theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) + facet_wrap(~ region, ncol = 1)+ labs(x = "NDVI", y = "dengue cases per 100,000")
p2 <- ggplot(tmp7, aes(x=NDVI, y=mal_inc)) + geom_point(size = 0.5, alpha = 1/5) + theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) + facet_wrap(~ region, ncol = 1)+ labs(x = "NDVI", y = "malaria cases per 100,000")

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/inc_NDVI_reg.png", width=8,height=4,units="in",res=1200)
plot_comb <- plot_grid(p1, p2, align='h', labels=c('(a)', '(b)'), label_size = 12)
plot_comb
dev.off()

ggplot(tmp7, aes(x=NDVI, y=den_inc)) + geom_point(size = 1, alpha = 1/10) + theme(axis.text=element_text(size=10),axis.title=element_text(size=10)) + labs(x = "NDVI", y = "dengue cases per 100,000") 


