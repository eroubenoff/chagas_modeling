# Climate data downloaded from: https://psl.noaa.gov/data/gridded/data.ghcncams.html
# Reading raster help from: https://rpubs.com/boyerag/297592

library(tidyverse)
library(ncdf4)
library(raster)

nc_data <- nc_open('climate_data/air.mon.mean.nc')
# Save the print(nc) dump to a text file
{
  sink('climate_data/air.mon.mean.metadata.txt')
  print(nc_data)
  sink()
}

lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")

head(lon) # look at the first few entries in the longitude vector

airtemp.array <- ncvar_get(nc_data, "air") # store the data in a 3-dimensional array
dim(airtemp.array) 
fillvalue <- ncatt_get(nc_data, "air", "missing_value")
fillvalue
nc_close(nc_data) 

# airtemp.array[airtemp.array == fillvalue$value] <- NA
airtemp.slice <- airtemp.array[, , 600] 
r <- raster(t(airtemp.slice), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
# r <- flip(r, direction='y')
plot(r)
