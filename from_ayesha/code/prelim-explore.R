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
#### load data


# area by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")

# pop by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")

#############################################################################
##### population density versus cases

# 2001
# dengue data
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_2001.Rdata")
deng <- data

# malaria data
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_2001.Rdata")
mal <- data

# urban pop by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned_1991-2000.Rdata")


pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/cases_vs_dem2001.pdf", width = 9, height = 6)
par(family = "serif")
dis <- c("deng", "mal")
for (d in 1:2) {
  data <- eval(parse(text = dis[d]))
  tmp <- data[,c("code", names(data)[which(grepl("Total", names(data)) == TRUE)])]
  names(tmp) <- c("code", "total")
  tmp2 <- left_join(tmp, area, by = c("code"))
  
  pop.tmp <- pop[,c("code", "2000")]
  names(pop.tmp) <- c("code", "pop")
  tmp3 <- left_join(tmp2, pop.tmp, by = c("code"))
  tmp3$pop <- as.numeric(as.character((tmp3$pop)))
  tmp3$area_km2 <- as.numeric(as.character((tmp3$area_km2)))
  tmp3$density <- tmp3$pop / tmp3$area_km2
  
  urban2 <- urban[,c("code", "urban2000")]
  tmp4 <- left_join(tmp3, urban2, by = c("code"))
  tmp4$urbanPop <- as.numeric(as.character(tmp4$urban2000))
  tmp4$urbanPerc <- (tmp4$urbanPop/tmp4$pop) * 100
  tmp4$urbanPerc[tmp4$urbanPerc > 100] <- 100
  
  plot.tmp <- as.data.frame(cbind(tmp4$code, tmp4$pop, tmp4$density, tmp4$urbanPerc, tmp4$total))
  names(plot.tmp) <- c("code","pop", "density","urban", "total")
  plot.tmp$total[is.na(plot.tmp$total)] <- 0
  if(d == 1) par(mfrow = c(2,3))
  plot(plot.tmp$density, plot.tmp$total, main = if(d == 1) "Dengue" else "Malaria", ylab = "cases",
       xlab = "density", pch = 20)
  plot(plot.tmp$pop, plot.tmp$total, xlab = "pop", ylab = "", pch = 20)
  plot(plot.tmp$urban, plot.tmp$total, xlab = "% urban", ylab = "", pch = 20)
  
}
dev.off()


#####

# 2010
# dengue data
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_2010.Rdata")
deng <- data

# malaria data
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_2010.Rdata")
mal <- data


# urbanizaton % by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned.Rdata")


pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/cases_vs_dem2010.pdf", width = 9, height = 6)
par(family = "serif")
dis <- c("deng", "mal")
for (d in 1:2) {
  data <- eval(parse(text = dis[d]))
  tmp <- data[,c("code", names(data)[which(grepl("Total", names(data)) == TRUE)])]
  names(tmp) <- c("code", "total")
  tmp2 <- left_join(tmp, area, by = c("code"))
  
  pop.tmp <- pop[,c("code", "2010")]
  names(pop.tmp) <- c("code", "pop")
  tmp3 <- left_join(tmp2, pop.tmp, by = c("code"))
  tmp3$pop <- as.numeric(as.character((tmp3$pop)))
  tmp3$area_km2 <- as.numeric(as.character((tmp3$area_km2)))
  tmp3$density <- tmp3$pop / tmp3$area_km2
  
  tmp4 <- left_join(tmp3, urban, by = c("code"))
  tmp4$Urban <- as.numeric(as.character(tmp4$Urban))
  
  plot.tmp <- as.data.frame(cbind(tmp4$code, tmp4$pop, tmp4$density, tmp4$Urban, tmp4$total))
  names(plot.tmp) <- c("code","pop", "density","urban", "total")
  plot.tmp$total[is.na(plot.tmp$total)] <- 0
  if(d == 1) par(mfrow = c(2,3))
  plot(plot.tmp$density, plot.tmp$total, main = if(d == 1) "Dengue" else "Malaria", ylab = "cases",
       xlab = "density", pch = 20)
  plot(plot.tmp$pop, plot.tmp$total, xlab = "pop", ylab = "", pch = 20)
  plot(plot.tmp$urban, plot.tmp$total, xlab = "% urban", ylab = "", pch = 20)
  
}
dev.off()


#############################################################################
# plot map of cases for each year

#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
municBR@data$order <- 1:nrow(municBR@data)


#### dengue
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,13,16)))
for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c("code", names(data)[which(grepl("Total", names(data)) == TRUE)])]
  names(tmp) <- c("code", "total")
  
  map.tmp <- municBR
  map.tmp@data <- merge(map.tmp@data, tmp,
                        by='code', all.x=T)
  
  myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd'))
  breaks = c(0,1,10,20,50,100,500,1000,3000, 10000, 100000000)
  
  
  for (i in 1:(length(myPalette) - 1)){
    map.tmp@data$color[(as.numeric(as.character(map.tmp@data$total))) >= breaks[i] & (as.numeric(as.character(map.tmp@data$total))) < breaks[i+1]] <- myPalette[i]
  }
  
  levs <- rep(NA, length(myPalette))
  for (i in 2:(length(myPalette) - 1)){
    levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
  }
  levs[10] <- "10000+"
  levs[1] <- "0"
  
  map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]
  
  
  if (y == 1) pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/dengSINAN.pdf", width = 7, height = 6)
  par(mar = c(0,0,1,0)); par(mfrow = c(1,1))
  plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, main = paste(years[y]))
  legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
  if (y == length(years)) dev.off()
}



#### malaria
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,14,17)))
for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c("code", names(data)[which(grepl("Total", names(data)) == TRUE)])]
  names(tmp) <- c("code", "total")
  
  map.tmp <- municBR
  map.tmp@data <- merge(map.tmp@data, tmp,
                        by='code', all.x=T)
  
  myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd')[4:9])
  breaks = c(0,1,5,10,20,40,80,100) 
  
  
  for (i in 1:(length(myPalette) - 1)){
    map.tmp@data$color[(as.numeric(as.character(map.tmp@data$total))) >= breaks[i] & (as.numeric(as.character(map.tmp@data$total))) < breaks[i+1]] <- myPalette[i]
  }
  
  levs <- rep(NA, length(myPalette))
  for (i in 2:(length(myPalette) - 1)){
    levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
  }
  levs[length(levs)] <- "100+"
  levs[1] <- "0"
  
  map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]
  
  if (y == 1) pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/malSINAN.pdf", width = 7, height = 6)
  par(mar = c(0,0,1,0)); par(mfrow = c(1,1))
  plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, main = paste(years[y]))
  legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
  if (y == length(years)) dev.off()
}



#############################################################################
# heatmap
library(fields)
library(gplots)

#### dengue
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,13,16)))
time <- seq(2001, 2013, 1/12)[1:144]
mtrx <- matrix(NA, 144, 27)
time.range <- seq(1,144, by = 12)

# load region masterfile for ordering municipalities
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")

# now fill in the matrix

for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c(2,4:15)]
  tmp2 <- inner_join(region_master, tmp, by = "code")
  tmp2[is.na(tmp2)]<- 0
  df = tmp2[,c(1,6:17)] %>% group_by(stateID) %>% summarise_each(funs(sum))
  
  mtrx[time.range[y]:(time.range[y]+11),1:27] <- t(as.matrix(df[,2:13]))
}
  
row.names(mtrx) <- time
col.labels <- unique(region_master$stateName)
row.labels <- rep("", length(time))
row.labels[seq(1,144, by = 12)] <- time[seq(1,144, by = 12)]

myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd'))
breaks = c(0,1,500,1000,5000, 10000,20000, 40000,60000,80000,100000)


dev.off()
pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/dengue_timeseries_states.pdf", width = 7, height = 6)

par(mfrow = c(1,1))
par(family="serif")
mar.default <- c(5.1, 4.1, 4.1, 2.1)
par(mar=c(5.1, 8, 4.1, 2.1))
cex = 1
image.plot(x=1:144,y=1:27, z=mtrx,xaxt="n",col=myPalette, breaks= breaks,
           cex=cex, cex.lab = cex, cex.axis = cex,cex.main=cex, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(1,1:27,labels=col.labels, xpd=TRUE, pos = 2,cex = cex)
axis(1, at=seq(1,144, by = 12),labels=time[seq(1,144, by = 12)],tick=F, cex = cex, cex.lab = cex, cex.axis = cex,cex.main=cex)
dev.off()



#### malaria
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,14,17)))
time <- seq(2001, 2015, 1/12)[1:168]
mtrx <- matrix(NA, 168, 27)
time.range <- seq(1,168, by = 12)

# load region masterfile for ordering municipalities
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")

# now fill in the matrix

for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c(2,4:15)]
  tmp2 <- inner_join(region_master, tmp, by = "code")
  tmp2[is.na(tmp2)]<- 0
  df = tmp2[,c(1,6:17)] %>% group_by(stateID) %>% summarise_each(funs(sum))
  
  mtrx[time.range[y]:(time.range[y]+11),1:27] <- t(as.matrix(df[,2:13]))
}

row.names(mtrx) <- time
col.labels <- unique(region_master$stateName)
row.labels <- rep("", length(time))
row.labels[seq(1,168, by = 12)] <- time[seq(1,168, by = 12)]

range(mtrx)
myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd'))
breaks = c(0,1,5,10,20,30,40,50,60,80,100)


dev.off()
pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/mal_timeseries_states.pdf", width = 7, height = 6)

par(mfrow = c(1,1))
par(family="serif")
mar.default <- c(5.1, 4.1, 4.1, 2.1)
par(mar=c(5.1, 8, 4.1, 2.1))
cex = 1
image.plot(x=1:168,y=1:27, z=mtrx,xaxt="n",col=myPalette, breaks= breaks,
           cex=cex, cex.lab = cex, cex.axis = cex,cex.main=cex, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(1,1:27,labels=col.labels, xpd=TRUE, pos = 2,cex = cex)
axis(1, at=seq(1,168, by = 12),labels=time[seq(1,168, by = 12)],tick=F, cex = cex, cex.lab = cex, cex.axis = cex,cex.main=cex)

dev.off()



#############################################################################
### heatmap ordered by places that had the biggest demographic changes
# 2000 to 2010

# area by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")

# pop by munic
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")



pop.tmp <- pop[,c("code", "2000", "2010")]
names(pop.tmp) <- c("code", "pop2000", "pop2010")
tmp <- left_join(area, pop.tmp, by = c("code"))
tmp$pop2000 <- as.numeric(as.character((tmp$pop2000)))
tmp$pop2010 <- as.numeric(as.character((tmp$pop2010)))

tmp$area_km2 <- as.numeric(as.character((tmp$area_km2)))
tmp$dens2000 <- tmp$pop2000 / tmp$area_km2
tmp$dens2010 <- tmp$pop2010 / tmp$area_km2

# urban pop by munic 2000
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned_1991-2000.Rdata")

urban2 <- urban[,c("code", "urban2000")]
tmp2 <- left_join(tmp, urban2, by = c("code"))
tmp2$urbanPop <- as.numeric(as.character(tmp2$urban2000))
tmp2$urban2000 <- (tmp2$urbanPop/tmp2$pop2000) * 100
tmp2$urban2000[tmp2$urban2000 > 100] <- 100


# urbanizaton % by munic 2010
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned.Rdata")

tmp3 <- left_join(tmp2, urban, by = c("code"))
tmp3$urban2010 <- as.numeric(as.character(tmp3$Urban))

diff.dat <- as.data.frame(cbind(tmp3$code, tmp3$pop2000, tmp3$pop2010.x, 
                                tmp3$dens2000, tmp3$dens2010, tmp3$urban2000, tmp3$urban2010))

names(diff.dat) <- c("code","pop2000","pop2010", "dens2000","dens2010","urban2000","urban2010")
diff.dat$urbanDiff <- diff.dat$urban2010 - diff.dat$urban2000
diff.dat$densDiff <- diff.dat$dens2010 - diff.dat$dens2000
diff.dat$popDiff <- diff.dat$pop2010 - diff.dat$pop2000


### map
#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
municBR@data$order <- 1:nrow(municBR@data)


map.tmp <- municBR
map.tmp@data <- merge(map.tmp@data, diff.dat,
                      by='code', all.x=T)
# urbanization

summary(map.tmp@data$urbanDiff, na.rm = T)
myPalette <- c(brewer.pal(11,'RdYlGn'))
breaks = c(-80,-60,-40,-20,-10,0,5,10,20,40,60,80)


for (i in 1:(length(myPalette) - 1)){
  map.tmp@data$color[(as.numeric(as.character(map.tmp@data$urbanDiff))) >= breaks[i] & (as.numeric(as.character(map.tmp@data$urbanDiff))) < breaks[i+1]] <- myPalette[i]
}

levs <- rep(NA, length(myPalette))
for (i in 2:(length(myPalette) - 1)){
  levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
}
levs[length(levs)] <- "80+"
levs[1] <- "less than -80"

map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]


pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/urbanDiff_2000_2010.pdf", width = 7, height = 6)
par(mar = c(0,0,1,0)); par(mfrow = c(1,1))
plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, 
     main = "% point change in urbanization")
legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
dev.off()

#density
summary(map.tmp@data$densDiff, na.rm = T)
myPalette <- c(brewer.pal(11,'RdYlGn'))
breaks = c(-300,-50,-20,-10,-5,0,5,10,20,100,500,2500)



for (i in 1:(length(myPalette) - 1)){
  map.tmp@data$color[map.tmp@data$densDiff >= breaks[i] & map.tmp@data$densDiff < breaks[i+1]] <- myPalette[i]
}

levs <- rep(NA, length(myPalette))
for (i in 2:(length(myPalette) - 1)){
  levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
}
levs[length(levs)] <- "2500+"
levs[1] <- "less than -300"

map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]


pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/densDiff_2000_2010.pdf", width = 7, height = 6)
par(mar = c(0,0,1,0)); par(mfrow = c(1,1))
plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, 
     main = "change in density (people / km2)")
legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
dev.off()

## pop size
summary(map.tmp@data$popDiff, na.rm = T)
myPalette <- c(brewer.pal(11,'RdYlGn'))
#breaks = c(seq(-2e+05, -1e+05, length.out = 5), 0, seq(1e+05, 1e+06, length.out = 6))

breaks = c(seq(-2e+05, -1e+05, length.out = 5), 0, 1000, 3000, 5000, 1e+05, 5e+05,1e+06)




for (i in 1:(length(myPalette) - 1)){
  map.tmp@data$color[map.tmp@data$popDiff >= breaks[i] & map.tmp@data$popDiff < breaks[i+1]] <- myPalette[i]
}

levs <- rep(NA, length(myPalette))
for (i in 2:(length(myPalette) - 1)){
  levs[i] <- paste0(breaks[i], " to ", breaks[i+1])
}
levs[length(levs)] <- "+1000000"
levs[1] <- "less than -200000"

map.tmp@data <- map.tmp@data[order(map.tmp@data$order),]


pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/popDiff_2000_2010.pdf", width = 7, height = 6)
par(mar = c(0,0,1,0)); par(mfrow = c(1,1))
plot(map.tmp, col=map.tmp@data$color, lty = 1, border = "grey", lwd = 0.3, 
     main = "change in population size")
legend("bottomright", fill = myPalette, legend = levs, col = myPalette, cex=0.8, bty = "n")
dev.off()



#### heatmap
topUrb <- diff.dat$code[order(-diff.dat$urbanDiff)][1:30]
topDens <- diff.dat$code[order(-diff.dat$densDiff)][1:30]
topPop <- diff.dat$code[order(-diff.dat$popDiff)][1:30]


# load region masterfile and ordering municipalities
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")
orderUrb <- region_master[region_master$code %in% topUrb,]
orderUrb <- orderUrb[match(topUrb, orderUrb$code),]

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,13,16)))
time <- seq(2001, 2013, 1/12)[1:144]
mtrx <- matrix(NA, 144, 30)
time.range <- seq(1,144, by = 12)


# now fill in the matrix

for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c(2,4:15)]
  tmp2 <- inner_join(orderUrb, tmp, by = "code")
  tmp2[is.na(tmp2)]<- 0
  df = tmp2[,c(6:17)] 
  
  mtrx[time.range[y]:(time.range[y]+11),1:30] <- t(as.matrix(df))
}

row.names(mtrx) <- time
col.labels <- orderUrb$code
row.labels <- rep("", length(time))
row.labels[seq(1,144, by = 12)] <- time[seq(1,144, by = 12)]

summary(mtrx)

myPalette <- c("#FFFFFF", brewer.pal(9,'YlOrRd'))
breaks = c(0,200,500,1000,5000, 10000,20000, 40000,60000,80000,100000)

breaks = c(0,5,10,20,50, 100,120, 140,160,180,200)


par(mfrow = c(1,1))
par(family="serif")
mar.default <- c(5.1, 4.1, 4.1, 2.1)
par(mar=c(5.1, 8, 4.1, 2.1))
cex = 1
image.plot(x=1:144,y=1:30, z=mtrx,xaxt="n",col=myPalette, breaks= breaks,
           cex=cex, cex.lab = cex, cex.axis = cex,cex.main=cex, xaxt = "n", yaxt = "n", xlab = "", ylab = "")
text(1,1:30,labels=col.labels, xpd=TRUE, pos = 2,cex = cex)
axis(1, at=seq(1,144, by = 12),labels=time[seq(1,144, by = 12)],tick=F, cex = cex, cex.lab = cex, cex.axis = cex,cex.main=cex)

#############################################################################
# urbanization versus dengue incidence

head(diff.dat)
tmp1 <- diff.dat
#dengue incidence in three year period (2001-2003 versus 2010-2012)
years <- c(2001:2003, 2010:2012)
for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c("code", names(data)[which(grepl("Total", names(data)) == TRUE)])]
  names(tmp) <- c("code", as.character(years[y]))
  tmp1 <- left_join(tmp1, tmp, by = "code")
  
}

tmp1$inc2001 <- ((tmp1$'2001' + tmp1$'2002' + tmp1$'2003') / (3 * tmp1$pop2000)) * 1000
tmp1$inc2010 <- ((tmp1$'2010' + tmp1$'2011' + tmp1$'2012') / (3 * tmp1$pop2010)) * 1000
tmp1$incDiff <- tmp1$inc2010 - tmp1$inc2001

summary(lm(tmp1$incDiff ~ tmp1$urbanDiff))
## NOT significant
#############################################################################


