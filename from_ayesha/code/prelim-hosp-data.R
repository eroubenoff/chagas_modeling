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

rm(list = ls())

#############################################################################
## Dengue
#############################################################################
# load data and get year and month from column names
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-PlaceOfHosp/DengueHOSP_Jan2008-Out2017.Rdata")
dates <- colnames(data)[4:(ncol(data)-1)]
months <- substr(dates, 6, 8)
years <- substr(dates, 1, 4)


#############################################################################
# plot map of cases for each year
numyears <- length(unique(years))
uniq.yr <- unique(years)


#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
municBR@data$order <- 1:nrow(municBR@data)


for (y in 1:numyears) {
  keep.cols <- which(substr(dates, 1, 4) == uniq.yr[y])
  tmp <- data[c(2,(3+keep.cols[1]):(3+keep.cols[length(keep.cols)]))]
  tmp$sum = rowSums(tmp[2:ncol(tmp)], na.rm = TRUE)
  tmp <- tmp[,c("code", "sum")]
  
  map.tmp <- municBR
  map.tmp@data <- merge(map.tmp@data, tmp,
                        by='code', all.x=T)
  
  myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd'))
  breaks = c(0,1,10,20,50,100,500,1000,3000, 10000, 100000000)
  
  
  for (i in 1:(length(myPalette) - 1)){
    map.tmp@data$color[(as.numeric(as.character(map.tmp@data$sum))) >= breaks[i] & (as.numeric(as.character(map.tmp@data$sum))) < breaks[i+1]] <- myPalette[i]
  }
  
  levs <- rep(NA, length(myPalette))
  for (i in 2:(length(myPalette) - 1)){
    levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
  }
  levs[10] <- "10000+"
  levs[1] <- "0"
  
  map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]
  
  
  if (y == 1) pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/dengHosp.pdf", width = 7, height = 6)
  par(mar = c(0,0,1,0))
  plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, main = paste(uniq.yr[y]))
  legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
  if (y == numyears) dev.off()
  


}




#############################################################################
## Malaria
#############################################################################
rm(list = ls())
# load data and get year and month from column names
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-PlaceOfHosp/MalariaHOSP_Jan2008-Out2017.Rdata")
dates <- colnames(data)[4:(ncol(data)-1)]
months <- substr(dates, 6, 8)
years <- substr(dates, 1, 4)


#############################################################################
# plot map of cases for each year
numyears <- length(unique(years))
uniq.yr <- unique(years)


#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
municBR@data$order <- 1:nrow(municBR@data)


for (y in 1:numyears) {
  keep.cols <- which(substr(dates, 1, 4) == uniq.yr[y])
  tmp <- data[c(2,(3+keep.cols[1]):(3+keep.cols[length(keep.cols)]))]
  tmp$sum = rowSums(tmp[2:ncol(tmp)], na.rm = TRUE)
  tmp <- tmp[,c("code", "sum")]
  
  map.tmp <- municBR
  map.tmp@data <- merge(map.tmp@data, tmp,
                        by='code', all.x=T)
  
  
  myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd')[4:9])
  breaks = c(0,1,5,10,20,40,80,100) 
  
  for (i in 1:(length(myPalette) - 1)){
    map.tmp@data$color[(as.numeric(as.character(map.tmp@data$sum))) >= breaks[i] & (as.numeric(as.character(map.tmp@data$sum))) < breaks[i+1]] <- myPalette[i]
  }
  
  levs <- rep(NA, length(myPalette))
  for (i in 2:(length(myPalette) - 1)){
    levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
  }
  levs[length(levs)] <- "100+"
  levs[1] <- "0"
  
  map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]
  
  
  if (y == 1) pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/malariaHosp.pdf", width = 7, height = 6)
  par(mar = c(0,0,1,0))
  plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, main = paste(uniq.yr[y]))
  legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
  if (y == numyears) dev.off()
  
  
  
}
#############################################################################
#### Plot time series
#############################################################################
# load data and get year and month from column names
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-PlaceOfHosp/DengueHOSP_Jan2008-Out2017.Rdata")
#dates <- colnames(data)[4:(ncol(data)-1)]
#months <- substr(dates, 6, 8)
#years <- substr(dates, 1, 4)

keycol <- "date"
valuecol <- "den"
gathercols <- colnames(data)[4:127]

hosp.long <- gather_(data, keycol, valuecol, gathercols) %>% dplyr::select(code, date, den) %>%
  mutate(year = as.numeric(substr(date,1,4)), 
         month = case_when(substr(date,6,8) == "Jan" ~ 1,
                           substr(date,6,8) == "Fev" ~ 2,
                           substr(date,6,8) == "Mar" ~ 3,
                           substr(date,6,8) == "Abr" ~ 4,
                           substr(date,6,8) == "Mai" ~ 5,
                           substr(date,6,8) == "Jun" ~ 6,
                           substr(date,6,8) == "Jul" ~ 7,
                           substr(date,6,8) == "Ago" ~ 8,
                           substr(date,6,8) == "Set" ~ 9,
                           substr(date,6,8) == "Out" ~ 10,
                           substr(date,6,8) == "Nov" ~ 11,
                           substr(date,6,8) == "Dez" ~ 12)) %>%
  dplyr::select(code, den, month, year)

head(hosp.long)

###### 1992 - 2007
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-PlaceOfHosp/DengueHOSP_1992-2007.Rdata")
head(data)
keycol <- "date"
valuecol <- "den"
gathercols <- colnames(data)[4:117]

hosp.long2 <- gather_(data, keycol, valuecol, gathercols) %>% dplyr::select(code, date, den) %>%
  mutate(year = as.numeric(substr(date,1,4)), 
         month = case_when(substr(date,6,8) == "Jan" ~ 1,
                           substr(date,6,8) == "Fev" ~ 2,
                           substr(date,6,8) == "Mar" ~ 3,
                           substr(date,6,8) == "Abr" ~ 4,
                           substr(date,6,8) == "Mai" ~ 5,
                           substr(date,6,8) == "Jun" ~ 6,
                           substr(date,6,8) == "Jul" ~ 7,
                           substr(date,6,8) == "Ago" ~ 8,
                           substr(date,6,8) == "Set" ~ 9,
                           substr(date,6,8) == "Out" ~ 10,
                           substr(date,6,8) == "Nov" ~ 11,
                           substr(date,6,8) == "Dez" ~ 12)) %>%
  dplyr::select(code, den, month, year) 

head(hosp.long2)

hosp.long.all <- rbind(hosp.long, hosp.long2)

# now merge all files
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/clim_all.Rdata")

hosp.long.all$year <- as.numeric(hosp.long.all$year)
tmp1 <- left_join(hosp.long.all, all_var, by = c("code", "year", "month"))

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
tmp6 <- tmp5 %>% group_by(code) %>% mutate(lag1_den=lag(den), lag1_pre = lag(pre), lag1_tmp = lag(tmp))

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

# create exposure variable
r = tmp7 %>% group_by(month) %>% summarise(r = sum(den, na.rm = TRUE)/sum(ps, na.rm = TRUE))
tmp7 <- left_join(tmp7, r, by = "month")
tmp7$exp.count <- tmp7$ps * tmp7$r

save(tmp7, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic_HOSPITAL.Rdata")


# plot cases
load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic_HOSPITAL.Rdata")
ts <- tmp7 %>% group_by(region, year, month) %>% 
  summarize(den = sum(den,na.rm = T), ps = sum(ps, na.rm = T)) %>%
  mutate(den_inc = den*100000/ps) %>% ungroup() %>%
  arrange(region, year, month) %>%
  mutate(time = rep(seq(1998, 2018, by = 1/12)[1:(length(seq(1998, 2018, by = 1/12))-3)], 5))





pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/timeseries_HOSPITAL.pdf", width=8,height=6)
ggplot(ts[ts$time <= 2015,], aes(x = time, y = den_inc, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
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



ggplot(ts, aes(x = time, y = den/1000, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_classic() +
  labs(x = "Month", y = "Monthly dengue cases (1000)",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())

ts2 <- tmp7 %>% group_by(region, year, month) %>% summarize(tmp = mean(tmp,na.rm = T), 
                                                     urban.perc = mean(urban.perc, na.rm = T),
                                                     gdp = mean(gdp, na.rm = T),
                                                     agric = mean(agric, na.rm = T),
                                                     NDVI = mean(NDVI, na.rm = T)) %>%
  arrange(region, year, month) %>% ungroup() %>%
  mutate(time =  rep(seq(1998, 2018, by = 1/12)[1:(length(seq(1998, 2018, by = 1/12))-3)], 5))

ggplot(ts2, aes(x = time, y = tmp, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_classic() +
  labs(x = "Month", y = "",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())

