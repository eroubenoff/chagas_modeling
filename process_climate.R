# Process climate data
# Data are downloaded from the NOAA PSL CPC
# https://psl.noaa.gov/data/gridded/data.cpc.globaltemp.html
# https://downloads.psl.noaa.gov/Datasets/cpc_global_temp/ 
# https://psl.noaa.gov/data/gridded/data.gpcc.html
# https://downloads.psl.noaa.gov/Datasets/cpc_global_precip/

library(tidyverse)
library(sf)
library(stars)
library(tmap)
# First read in Brazil shapefile
setwd("~/chagas_modeling")
load("chagas_data.Rdata")

# Get bounding box for Brazil

# Loop through each file and subset it for Brazil's boundaries
# Overwrite the original object-- cuts down on memory
path <- "climate_data/CPC"
flist <- list.files(path = path, pattern = ".nc")
flist
newpath <- "climate_data/CPC/Brazil"

process_climate <- function(f, br_shp) {
  br_bb <- st_bbox(br_shp)
  r <- read_ncdf(f, crs = st_crs(br_shp))
  r <- st_crop(r, br_bb)
  r_poly <- st_as_sf(r, as_points=FALSE)
  r_poly <- units::drop_units(r_poly) 
  
  r_new <- br_shp %>% 
    select(muni_code) %>%
    st_join(r_poly) %>% 
    mutate(area = units::drop_units(st_area(.))) %>%
    st_drop_geometry() %>%
    group_by(muni_code) %>%
    summarize(across(everything(), ~weighted.mean(., area)))
  
  return(r_new)
}




overwrite <- FALSE
for (f in flist) {
  out_path <- str_replace_all(f, "\\.", "_")
  out_path <- str_replace(out_path, "_nc", ".RDS")
  out_path <- file.path(newpath, out_path)
  
  if (!overwrite) {
    if (file.exists(file.path(newpath, out_path))) {
      message("Skipping ", f)
      next
    }
  }
  
  var <- process_climate(file.path(path,f), br_shp)
  saveRDS(var, file=out_path)
}


