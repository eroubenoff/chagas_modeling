require(raster)
require(sp)
require(MODIS)
library(rgdal)
library(gdalUtils)
library(raster)
library(rgeos)
library(velox)
library(rts)

################## DIDN'T WORK
## download files
years <- 2013:2017
months <- c(paste0("0",1:9), paste0(10:12))
for (y in 1:length(years)) {
  for (m in 1:length(months)) {
    page <- readLines(paste0("https://e4ftl01.cr.usgs.gov/MOLT/MOD13C2.006/", years[y],".", months[m],".01/"))
    setwd(paste0("~/Google Drive/Dengue/Dengue-Malaria/MODIS/", years[y],"/", sep = ""))
    filename <- substr(page[33],93,130)
    url <- paste0("https://e4ftl01.cr.usgs.gov/MOLT/MOD13C2.006/", years[y],".", months[m],".01/",filename)
    downloadFile(url, filename, username = "amahmud", password ="Full$ail4", binary = TRUE,mode = "curl", quiet = FALSE, cacheOK = T)
    }
}

gdalinfo(filename)
url 
##################
rm(list = ls())

# get municipality shapefile
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/br_municipios/")
munic <- readOGR(dsn=".", layer="BRMUE250GC_SIR")

year <- "2016"
setwd(paste("~/Google Drive/Dengue/Dengue-Malaria/MODIS/", year, sep = ""))

hdf_path <- paste("/Users/mahmud/Google Drive/Dengue/Dengue-Malaria/MODIS/", year, sep = "")
hdf_files <- list.files(hdf_path, pattern="^MOD13C2.*\\.hdf$")
stopifnot(length(hdf_files) > 0)

hdf_dir <- list.files(hdf_path, pattern="hdf$", full.names=FALSE) #create list of names with hdf extension
#Create full directory name of files
fhdf <-list.files(hdf_path, pattern="hdf$", full.names=TRUE)

for (i in 1:length(hdf_files)) {
  infile <- hdf_files[i]
  
  namefi <- substring(hdf_dir[i],1,nchar(hdf_dir[i])-4)
  outfile <- paste0(hdf_path,"/", namefi,".tif")
  gdal_translate(fhdf[i], outfile, sd_index=2)
  #index = 1 is NDVI, index = 2 is EVI
  #The Enhanced Vegetation Index (EVI) improves on the NDVI. 
  
  tifFile <- substr(outfile, 61,nchar(outfile))

  vi <- preStack(tifFile)
  s <- stack(vi)
  
  s2 <- s*0.0001 # Rescale the downloaded Files with the scaling factor
  # see here: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13c2_v006
  
  sp <- spTransform(munic,CRS(proj4string(s)))
  
  v1 <- velox(s2*0.0001, extent = c(-180, 180, -90, 90), res = c(0.05, 0.05), crs = "+proj=longlat +ellps=clrk66 +no_defs")
  ext <- v1$extract(sp = sp, fun = mean)
  sp@data$NDVI <- ext
  
  sp@data$code_long = as.character(unique(munic$CD_GEOCMU))
  sp@data$code = substr(sp@data$code_long, 1,6)
  
  data <- as.data.frame(sp@data)
  
  year = substr(tifFile, 10,13)
  month = i
  #ceiling(as.numeric(substr(tifFile, 14,16))/30) #only works for leap years
  #Julian day calendar: https://landweb.modaps.eosdis.nasa.gov/browse/calendar.html
  
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_",year,"_",month,".Rdata", sep = "")
  save(data, file = filepath)
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_Spatial_NDVI",year,"_",month,".Rdata", sep = "")
  save(sp, file = filepath)
  
  print(i)
}



################## FOR 2012
rm(list = ls())

# get municipality shapefile
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/br_municipios/")
munic <- readOGR(dsn=".", layer="BRMUE250GC_SIR")

year <- "2012"
setwd(paste("~/Google Drive/Dengue/Dengue-Malaria/MODIS2/", year, sep = ""))

hdf_path <- paste("/Users/mahmud/Google Drive/Dengue/Dengue-Malaria/MODIS2/", year, sep = "")
hdf_files <- list.files(hdf_path, pattern="^MOD13C2.*\\.hdf$")
stopifnot(length(hdf_files) > 0)

hdf_dir <- list.files(hdf_path, pattern="hdf$", full.names=FALSE) #create list of names with hdf extension
#Create full directory name of files
fhdf <-list.files(hdf_path, pattern="hdf$", full.names=TRUE)

for (i in 1:length(hdf_files)) {
  infile <- hdf_files[i]
  
  namefi <- substring(hdf_dir[i],1,nchar(hdf_dir[i])-4)
  outfile <- paste0(hdf_path,"/", namefi,".tif")
  gdal_translate(fhdf[i], outfile, sd_index=2)
  #index = 1 is NDVI, index = 2 is EVI
  #The Enhanced Vegetation Index (EVI) improves on the NDVI. 
  
  tifFile <- substr(outfile, 62,nchar(outfile))
  
  vi <- preStack(tifFile)
  s <- stack(vi)
  
  s2 <- s*0.0001 # Rescale the downloaded Files with the scaling factor
  # see here: https://lpdaac.usgs.gov/dataset_discovery/modis/modis_products_table/mod13c2_v006
  
  sp <- spTransform(munic,CRS(proj4string(s)))
  
  v1 <- velox(s2*0.0001, extent = c(-180, 180, -90, 90), res = c(0.05, 0.05), crs = "+proj=longlat +ellps=clrk66 +no_defs")
  ext <- v1$extract(sp = sp, fun = mean)
  sp@data$NDVI <- ext
  
  sp@data$code_long = as.character(unique(munic$CD_GEOCMU))
  sp@data$code = substr(sp@data$code_long, 1,6)
  
  data <- as.data.frame(sp@data)
  
  year = substr(tifFile, 10,13)
  month = i
  #ceiling(as.numeric(substr(tifFile, 14,16))/30) #only works for leap years
  #Julian day calendar: https://landweb.modaps.eosdis.nasa.gov/browse/calendar.html
  
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_NDVI_",year,"_",month,".Rdata", sep = "")
  save(data, file = filepath)
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/modis/cleaned_Spatial_NDVI",year,"_",month,".Rdata", sep = "")
  save(sp, file = filepath)
  
  print(i)
}



