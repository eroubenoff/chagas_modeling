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

# Then, loop through all of the climate files
path <- "climate_data/CPC/Brazil"
flist <- list.files(path=path, full.names = TRUE)

# Read in each file for each variable and create a long data frame
precip_f <- flist %>% str_subset("precip")
tmin_f <- flist %>% str_subset("tmin")
tmax_f <- flist %>% str_subset("tmax")

precip_df <- vector("list", length(precip_f))
for (i in 1:length(precip_f)) {
  precip_df[[i]] <- read_rds(precip_f[i]) %>% st_drop_geometry()
}
precip_df <- reduce(precip_df, full_join, by="muni_code")
precip_df  <- precip_df %>%
  pivot_longer(-muni_code, names_to = "date", values_to = "precip")
precip_df <- precip_df %>% mutate(date = lubridate::as_date(date))


tmin_df <- vector("list", length(tmin_f))
for (i in 1:length(tmin_f)) {
  tmin_df[[i]] <- read_rds(tmin_f[i]) %>% st_drop_geometry()
}
tmin_df <- reduce(tmin_df, full_join, by="muni_code")
tmin_df  <- tmin_df %>%
  pivot_longer(-muni_code, names_to = "date", values_to = "tmin")
tmin_df <- tmin_df %>% mutate(date = lubridate::as_date(date))


tmax_df <- vector("list", length(tmax_f))
for (i in 1:length(tmax_f)) {
  tmax_df[[i]] <- read_rds(tmax_f[i]) %>% st_drop_geometry()
}
tmax_df <- reduce(tmax_df, full_join, by="muni_code")
tmax_df  <- tmax_df %>%
  pivot_longer(-muni_code, names_to = "date", values_to = "tmax")
tmax_df <- tmax_df %>% mutate(date = lubridate::as_date(date))


# Join them together
climate_df <- reduce(list(precip_df, tmin_df, tmax_df), left_join, by = c("muni_code", "date"))
rm(precip_df, tmin_df, tmax_df)
gc()


climate_df <- climate_df %>%
  mutate(date = lubridate::as_date(date))


# // Earth Engine Code: 
# var dataset = ee.ImageCollection('MODIS/061/MOD13A2')
# .filter(ee.Filter.date('2000-01-01', '2018-12-31'));
# var ndvi = dataset.select('NDVI');
# 
# 
# var ndviVis = {
#   min: 0.0,
#   max: 9000.0,
#   palette: [
#     'FFFFFF', 'CE7E45', 'DF923D', 'F1B555', 'FCD163', '99B718', '74A901',
#     '66A000', '529400', '3E8601', '207401', '056201', '004C00', '023B01',
#     '012E01', '011D01', '011301'
#   ],
# };
# Map.setCenter(6.746, 46.529, 2);
# // Map.addLayer(ndvi, ndviVis, 'NDVI');
# 
# var majority_value_per_poly = dataset.map(function(image) {
#   return image
#   .reduceRegions({
#     collection: br_shp,   // a shapefile of polygons ingested as asset
#     reducer: ee.Reducer.mean() //, // use mode to find the majority land cover per region (polygon)
#     // scale:100  // target pixel size
#   });
# 
# });
# 
# // print(majority_value_per_poly.flatten());
# // Export
# Export.table.toDrive({
#   collection: majority_value_per_poly.flatten(),  // <----- flatten()
#   description: 'br_ndvi',
#   fileFormat: 'CSV' ,
#   selectors: ['system:index', 'muni_code', 'NDVI', 'EVI']
# });








# Load the NDVI
NDVI <- read_csv("br_ndvi.csv")
NDVI <- NDVI %>% 
  separate_wider_delim(`system:index`, delim = "_", names = c("year", "month", "day", "idx"))

NDVI <- NDVI %>%
  mutate(date = lubridate::ymd(paste(year, month, day))) %>%
  select(muni_code, date, NDVI, EVI)

# NDVI is every 16 days but other climate vars aren't; 
# aggregate the other vars to the same time sequence as NDVI so can join

dr <- NDVI %>% pull(date) %>% unique() %>% sort()

climate_df <- climate_df %>%
  mutate(dr = cut(lubridate::as_date(date), breaks = dr, include.lowest = TRUE)) %>%
  group_by(dr) %>%
  summarize(max_tmax = max(tmax, na.rm=TRUE),
            min_tmax = min(tmax, na.rm=TRUE),
            mean_tmax = mean(tmax, na.rm=TRUE),
            tmax_gt_30 = sum(tmax > 30, na.rm=TRUE),
            tmax_gt_35 = sum(tmax > 35, na.rm=TRUE),
            max_tmin = max(tmin, na.rm=TRUE),
            min_tmin = min(tmin, na.rm=TRUE),
            mean_tmin = mean(tmin, na.rm=TRUE),
            max_precip = max(precip, na.rm=TRUE),
            min_precip = min(precip, na.rm=TRUE),
            mean_precip = mean(precip, na.rm=TRUE),
            precip_gt_
            )


