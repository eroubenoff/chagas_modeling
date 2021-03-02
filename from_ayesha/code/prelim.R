library(dplyr)

#Summary stats - dengue and malaria distribution overlap


load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN/DengueSINAN_2012.Rdata")
months <- colnames(data)[4:15]
months <- substr(months, 1, 3)
deng <- data

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN/MalariaSINAN_2012.Rdata")
mal <- data

mun.list <- list()
for (m in 1:length(months)) {
  mth <- paste(months[m], "_d",sep = "")
  tmp <- deng[c("code", mth)]
  
  mth <- paste(months[m], "_m",sep = "")
  tmp2 <-  mal[c("code", mth)]
  
  joined <- left_join(tmp, tmp2, by = "code")
  joined[is.na(joined)] <- 0
  
  joined$both <- ifelse(joined[,2] != 0 & joined[,3] != 0, 1, 0)
  sum(joined$both)
  mun.list[[m]] <- joined$code[which(joined$both == 1)]
  
}

#https://dioferrari.wordpress.com/2014/11/27/plotting-maps-using-r-example-with-brazilian-municipal-level-data/
########################################################
# map the states with overlap

require(maptools)
require(descr)
require(RColorBrewer)
require(plotrix)

#load shapefiles
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")

head(municBR@data)
municBR@data$code <- substr(municBR@data$codigo_ibg, 1,6)
length(unique(municBR@data$code))
length(unique(municBR@data$codigo_ibg))

for (i in 1:12) {
  tmp <- as.data.frame(cbind(mun.list[[i]], rep(1, length(mun.list[[i]]))))
  names(tmp) <- c("code", "present")
  
  map.tmp <- municBR
  map.tmp@data <- merge(map.tmp@data, tmp,
                        by='code', all.x=T)
  
  map.tmp@data$present[is.na(map.tmp@data$present)] <- 0
  map.tmp@data$colors <- ifelse(map.tmp@data$present == 1, "blue", "white")
  
  #table(municBR@data$present)
  
  if (i == 1) par(mfrow = c(4,3), family = "serif", mar = c(1,1,1,1))
  plot(map.tmp, col = map.tmp@data$colors,lty=0)
  plot(statesBR, add=T)
  #plot(regionsBR, add=T, lwd=3)
  
}



