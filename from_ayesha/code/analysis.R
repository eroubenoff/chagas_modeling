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

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_2012.Rdata")
months <- colnames(data)[4:15]
months <- substr(months, 1, 3)
deng <- data

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_2012.Rdata")
mal <- data

dat_long <- NA
for (m in 1:length(months)) {
  mth <- paste(months[m], "_d",sep = "")
  tmp <- deng[c("code", mth)]
  
  mth <- paste(months[m], "_m",sep = "")
  tmp2 <-  mal[c("code", mth)]
  
  joined <- left_join(tmp, tmp2, by = "code")
  joined[is.na(joined)] <- 0
  joined$month = m
  names(joined)[2:3] <- c("den", "mal")
  
  joined$both <- ifelse(joined[,2] != 0 & joined[,3] != 0, 1, 0)
  dat_long <- rbind(dat_long, joined)

}
dat_long <- dat_long[-1, ]

### Now bring in covariates
## climate data
## set variable name
var = c("pre", "tmn", "tmp", "tmx")

## set time range
range = c("2012-01-01", "2013-01-01")

rm(long, tmp)
for (v in 1:length(var)) {
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/climdata/cleaned_",var[v],"_",range[1],"_",range[2],".Rdata", sep = "")
  load(filepath)
  
  # reshape to long format
  keycol <- "date"
  valuecol <- var[v]
  gathercols <- paste0("Y2012M", seq(1:12))
  long <- gather_(data, keycol, valuecol, gathercols)
  long$year <- as.numeric(substr(long$date,2,5))
  long$month <- as.numeric(substr(long$date,7,8))
  long$code <- as.numeric(long$code)
  if (v == 1) tmp <- long
  if(v>= 2) {
    all_var <- left_join(long, tmp, by = c("code_long", "code", "date", "year", "month"))
    tmp <- all_var
  }
}


head(all_var)

filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_2012_",1,".Rdata", sep = "")
load(filepath)
tmp <- data
tmp$month <- 1
tmp$year <- 2012
## modis data
for (m in 2:12){
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_2012_",m,".Rdata", sep = "")
  load(filepath)
  data$month = m
  data$year = 2012
  tmp <- rbind(tmp, data)
  rm(data)
  
}
mod <- tmp
head(mod)
mod$code <- as.numeric(mod$code)

all_var <- left_join(all_var, mod, by = c("code", "month", "year"))

head(all_var)

# merge with data
dat <- left_join(dat_long, all_var, by = c("code", "month"))
dat$NDVI <- as.numeric(dat$NDVI)

# add in populaton, urban-rural and density
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned.Rdata")
dat_1 <- left_join(dat, urban, by = c("code"))
dat_1$Urban <- as.numeric(as.character((dat_1$Urban)))

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")
dat_2 <- left_join(dat_1, area, by = c("code"))
dat_2$area_km2 <- as.numeric(as.character((dat_2$area_km2)))


load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")
pop2012 <- pop[,c("code", "2012")]
names(pop2012) <- c("code", "pop2012")
dat_3 <- left_join(dat_2, pop2012, by = c("code"))
dat_3$density <- as.numeric(dat_3$pop2012) / dat_3$area_km2


fit <- glm(both ~ NDVI + tmx+ tmn + pre +  as.factor(month) + density + Urban, data = dat_3, family = "binomial")
summary(fit)
table(dat$both)

gam.fit <- gam(both ~ s(tmp) + s(pre) + as.factor(month), data = dat_3, family = "binomial")
summary(gam.fit)

dev.off()
plot(gam.fit, se = TRUE)

############################################################
#### map out where cases are 


#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")

head(municBR@data)
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
length(unique(municBR@data$code))
length(unique(municBR@data$codigo_ibg))

tmp <- dat_3[,c("code","den")]
head(tmp)
map.deng <- group_by(tmp, code) %>% summarise(sum = sum(den))

map.tmp <- municBR
map.tmp@data <- merge(map.tmp@data, map.deng,
                      by='code', all.x=T)
head(map.tmp@data)
summary(map.tmp@data$sum)


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

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/dengue_2012.pdf", width = 7, height = 6)
plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3)
legend("bottomleft", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
dev.off()

##use ggplot
# map.tmp@data$id = rownames(map.tmp@data)
# map.points <- fortify(map.tmp, region = "id")
# map.df = left_join(map.points, map.tmp@data, by="id")
# 
# 
# ggplot(map.df) + 
#   aes(long,lat,group=group,fill=sum) + 
#   geom_polygon() +
#   geom_path(color="white") +  scale_fill_gradientn(colours = myPaletteBlue,
#                                                    breaks = c(0,5,10,20,80, 150, 500, 1000, 5000),
#                                                    trans = "log10")
# 



#####
tmp <- dat_3[,c("code","mal")]
head(tmp)
map.mal <- group_by(tmp, code) %>% summarise(sum = sum(mal))

mal <- as.data.frame(map.mal[which(map.mal$sum > 0), ])

map.tmp <- municBR
map.tmp@data <- merge(map.tmp@data, map.mal,
                      by='code', all.x=T)
head(map.tmp@data)
summary(map.tmp@data$sum)

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

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/malaria_2012.pdf", width = 7, height = 6)
plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3)
legend("bottomleft", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
dev.off()

#####
tmp <- dat_3[,c("code","both")]
head(tmp)
map.mal <- group_by(tmp, code) %>% summarise(sum = sum(both))

map.tmp <- municBR
map.tmp@data <- merge(map.tmp@data, map.mal,
                      by='code', all.x=T)
head(map.tmp@data)
map.tmp@data$color[(as.numeric(as.character(map.tmp@data$sum))) > 0] <- "blue"


plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey")
