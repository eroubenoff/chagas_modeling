#-------------------------------------------------------------------------------
# 01_load_data.R
# This script loads the data and saves it as a R binary called 
#   "chagas_data.Rdata"
#-------------------------------------------------------------------------------


# Chagas Modeling 
library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)


#-------------------------------------------------------------------------------
#### Load Data ####
#-------------------------------------------------------------------------------

setwd("~/90days/eroubenoff/chagas")
load("./chagas_data/chagas_municipio.Rdata")
br_shp <- st_read("./from_ayesha/municipios_2010.shp", "municipios_2010")
load("./from_ayesha/cleaned_data/pop_all.Rdata")

#-------------------------------------------------------------------------------
#### Extract code from shapefile ####
#-------------------------------------------------------------------------------
br_shp <- br_shp %>%
  mutate(muni_code = as.numeric(substr(codigo_ibg, 1, 6)))





#-------------------------------------------------------------------------------
#### Rectangularize the Chagas data ####
#-------------------------------------------------------------------------------
# Pivot and extract muni code
chagas_municipio <- chagas_municipio %>%
  pivot_longer(cols = 2:14, names_to = "month", values_to = "count") %>%
  mutate(muni_code = str_sub(`Município de residência`, 1, 6),
         muni_code = as.numeric(muni_code))

# Create df
chagas_arr <- expand_grid(
  municipio = unique(br_shp$muni_code),
  year = unique(chagas_municipio$year),
  month = unique(chagas_municipio$month)
)

# Join df with data
chagas_arr <- left_join(chagas_arr,
                        chagas_municipio %>% select(
                          municipio = muni_code,
                          year, month, count),
)

# Add in UF 
chagas_arr <- left_join(chagas_arr,
                        br_shp %>%
                          st_drop_geometry() %>%
                          select(municipio = muni_code, uf))

# Set NAs to 0
chagas_arr[is.na(chagas_arr)] <- 0

# Lookup indices for municipios and UFs
chagas_arr <- left_join(
  chagas_arr,
  chagas_arr %>% 
    select(uf, municipio) %>%
    distinct() %>%
    arrange(uf, municipio) %>%
    group_by(uf) %>%
    mutate(uf_no = cur_group_id(),
           muni_no = 1:n())
) 



# Cast to array
# Need to order the month column like: factor(ColumnName, levels = unique(ColumnName)
chagas_count <- acast(chagas_arr, formula = 
                        uf_no~muni_no~year~factor(month, levels = unique(month)), 
                      value.var = "count",
                      fill = 0)
# dim(chagas_count)
# dimnames(chagas_count)

# Vector of the number of municipalities within each UF
N <- chagas_arr %>%
  select(uf_no, muni_no) %>%
  distinct() %>%
  group_by(uf_no) %>%
  slice_max(muni_no, n = 1) %>%
  pull(muni_no)





#-------------------------------------------------------------------------------
#### Rectangularize the Pop data ####
#-------------------------------------------------------------------------------
# Create pop df
chagas_pop <- chagas_arr %>%
  select(municipio, year, uf_no, muni_no) %>%
  distinct()

# Join df with data
chagas_pop <- chagas_pop %>%
  left_join(pop.all %>% select(municipio = code, year, ps))

# Cast matrix
chagas_pop <- acast(chagas_pop, formula = uf_no~muni_no~year, value.var = "ps", fill = 0)
# dim(chagas_pop)
# dimnames(chagas_pop)

# The population data only go through 2015 so for now we need to limit
# the incidence data to 2015 (NTS: get the rest of the pop data)
chagas_count <- chagas_count[, , 0:15, 1:12]
chagas_pop <- chagas_pop[, ,1:15] 

# See which data are missing by lining up the pop and the count data
if (FALSE) {
  for (u in 1:dim(chagas_pop)[1]) {
    message(u, "\n")
    for (m in 1:dim(chagas_pop)[2]) {
      for (y in 1:dim(chagas_pop)[3]) {
        if (is.na(chagas_pop[u, m, y]) & all(is.na(chagas_count[u, m, y, ]))){
          next
        }
        if (!is.na(chagas_pop[u, m, y]) & all(!is.na(chagas_count[u, m, y, ]))){
          next
        }
        else {
          message("Mismatch at: UF ", u, ", M ", m, ", Y ", y, "missmatch \n")
        }
      }
    }
  }
}

# There are a few missing populations
chagas_pop[12, 35, c(1,2,3)] <- rep(chagas_pop[12, 35, 4], 3)
chagas_pop[13, 51, c(1,2,3)] <- rep(chagas_pop[13, 51, 4], 3)
chagas_pop[13, 52, c(1,2,3)] <- rep(chagas_pop[13, 52, 4], 3)
chagas_pop[17, 14, c(1,2,3)] <- rep(chagas_pop[17, 14, 4], 3)






#-------------------------------------------------------------------------------
#### Write out object ####
#-------------------------------------------------------------------------------

save(
  chagas_count,
  chagas_pop,
  N,
  br_shp,
  chagas_arr,
  file = "chagas_data.Rdata"
)



