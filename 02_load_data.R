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
if (!require("devtools")) {
  install.packages("devtools")
}
# devtools::install_github("stan-dev/bayesplot")
library(bayesplot)
library(coda)


setwd("~/chagas_modeling")


#-------------------------------------------------------------------------------
#### Load Data ####
#-------------------------------------------------------------------------------

# setwd("~/90days/eroubenoff/chagas")
# load("./chagas_data/chagas_municipio.Rdata")
# load("./from_ayesha/cleaned_data/pop_all.Rdata")

#-------------------------------------------------------------------------------
#### Extract code from shapefile ####
#-------------------------------------------------------------------------------
# br_shp <- st_read("./from_ayesha/municipios_2010.shp", "municipios_2010")
br_shp <- geobr::read_municipality(year = 2015)
br_shp <- st_transform(br_shp, 4326)
# Drop 2 islands
to_drop <- c(1526, 3500)
br_shp <- br_shp[-to_drop, ]


## Load population
br_pop <- read_csv("br_pop.csv", skip = 3) %>% drop_na() # Last few lines are footnotes
br_pop <- br_pop %>% 
  rename(code_muni = 1, name = 2)  %>%
  mutate(code_muni = as.numeric(code_muni)) %>% 
  rename_with(~paste0("pop_", .), starts_with("2")) %>%
  mutate(across(starts_with("pop"), ~as.numeric(.)))

# This is missing the 2010 census, which is downloaded separately
flist <- list.files("br_2010_census_por_municipios", full.names = TRUE)
flist <- flist %>% str_subset("Tabela 1\\.")

br_2010 <- map_dfr(flist, ~{
  readxl::read_excel(., skip = 6) %>% 
    select(name = 1, pop = 2, code_muni = 11)
})

br_2010 <- br_2010 %>% filter(name != "Total")


# Join the 2010 census to the br_pop
br_pop <- br_pop %>% 
  left_join(br_2010 %>% select(-name, pop_2010 = pop, code_muni), by = c("code_muni") )


# Also missing 2007. Calculate 2007 as the average of 2006 and 2008.
br_pop <- br_pop %>%
  mutate(pop_2007 = round((pop_2006 + pop_2008)  / 2))


# Sort the years 
br_pop <- br_pop %>% select(code_muni, name, paste0("pop_", 2001:2021))



# There are a handful of municipalities missing data. Fill them in 
# with the closest year (I assume these were municipalities that 
# were created during this period)
br_pop[!complete.cases(br_pop),]
for (i in 21:3) {
  missing <- which(is.na(br_pop[,i]))
  br_pop[missing, i] <- br_pop[missing, i+1]
}
br_pop[!complete.cases(br_pop),]

# Check that all of these are in the shapefile and the CD dataset
which(!br_pop$code_muni %in% br_shp$code_muni)

# Join the other set of codes
br_pop <- full_join(br_shp, br_pop) 

# Two munis are missing population data: Lagoa Mirim and Lagoa Dos Patos. Drop
br_pop[!complete.cases(br_pop %>% st_drop_geometry()),]

br_pop <- br_pop[complete.cases(br_pop %>% st_drop_geometry()),]

br_pop


#-------------------------------------------------------------------------------
## Load Chagas counts
#-------------------------------------------------------------------------------
l <- list.files("chagas_data/from_tabnet", full.names = TRUE)
chagas_tabnet <- map_dfr(l, ~{
  df <- read_delim(., delim = ";", na = "-")
  year <- str_sub(., -8, -5)
  level <- if_else(str_detect(., "UF"), "UF", "municipio")
  df <- df %>%
    mutate(year = year,
           level = level) 
})

chagas_tabnet <- chagas_tabnet %>%
  filter(`Município de residência` != "Total")
chagas_municipio <- filter(chagas_tabnet, level == "municipio")

# Total cases: 
sum(chagas_municipio$Total)

# Take the annual summary:
chagas_municipio <- chagas_municipio %>%
  select('name' = `Município de residência`, Total, year) %>%
  mutate(muni_code = str_sub(name, 1, 6),
         muni_code = as.numeric(muni_code),
         name = str_sub(name, 8, -1))


# Pivot the table to wide format
chagas_municipio <- chagas_municipio %>% 
  pivot_wider(names_from = year, values_from = Total, values_fill = 0)

chagas_municipio <- chagas_municipio %>%
  rename_with(~paste0("chagas_", .), starts_with("2"))

nrow(chagas_municipio)

chagas_municipio



#-------------------------------------------------------------------------------
# Join Chagas data with the Pop data
#-------------------------------------------------------------------------------

# The population and shapefiles have a 7-digit code, but the Chagas 
# data has a 6-digit code. This thread explains it a little. 
# https://forum.ipums.org/t/municipalities-in-1991-brazilian-census/2847
# Basically can just disregard the 7th digit in the code

nDigits <- function(x) nchar( trunc( abs(x) ) )
nDigits <- Vectorize(nDigits)
nDigits(chagas_municipio$muni_code) %>% summary()
nDigits(br_pop$code_muni) %>% summary()

sum(br_pop$code_muni %in% chagas_municipio$muni_code)

br_pop <- br_pop %>% mutate(muni_code = as.numeric(str_sub(as.character(code_muni), 1, 6)))

sum(!chagas_municipio$muni_code %in% br_pop$muni_code)
which(!chagas_municipio$muni_code %in% br_pop$muni_code)
# Only one fails, so drop

br <- left_join(br_pop, chagas_municipio, by = c("muni_code"))

# Check that the names are close enough
# br %>% st_drop_geometry() %>% select(name.x, name.y) %>% View()

br <- br %>% select(-name.y) %>% rename(name = name.x)
br[is.na(br)] <- 0

#-------------------------------------------------------------------------------
#### Split the data into separate arrays ####
#-------------------------------------------------------------------------------

chagas_arr <- br %>%
  st_drop_geometry() %>%
  select(muni_code, contains("chagas")) %>% 
  pivot_longer(-muni_code) %>%
  mutate(year = as.numeric(str_sub(name, 8, -1))) %>%
  select(-name) %>%
  acast(year ~ muni_code, value.var = "value")

str(chagas_arr)
sum(chagas_arr)

pop_arr <- br %>%
  st_drop_geometry() %>%
  select(muni_code, contains("pop")) %>% 
  mutate(across(everything(), as.numeric)) %>%
  pivot_longer(-muni_code) %>%
  mutate(year = as.numeric(str_sub(name, 5, -1))) %>%
  filter(year %in% 2001:2019) %>%
  select(-name)  %>%
  acast(year ~ muni_code, value.var = "value")

str(pop_arr)
sum(pop_arr)

# Confirm that the dimensions are right
all(dimnames(chagas_arr)[[1]] == dimnames(pop_arr)[[1]])
all(dimnames(chagas_arr)[[2]] == dimnames(pop_arr)[[2]])


#-------------------------------------------------------------------------------
## Read in covariates 
#-------------------------------------------------------------------------------

covariate_df <- readRDS("bioclimate.RDS")

covariate_arr <- covariate_df %>% 
  mutate(year = lubridate::year(date)) %>% 
  filter(year >= 2001) %>% 
  filter(muni_code %in% pull(br,muni_code)) %>% 
  select(muni_code, year, contains("PC"))  %>%
  pivot_longer(-c(muni_code, year)) %>% 
  acast(year ~ muni_code ~ name)

all(dimnames(chagas_arr)[[1]] == dimnames(covariate_arr)[[1]][1:19])
all(dimnames(chagas_arr)[[2]] == dimnames(covariate_arr)[[2]])
covariate_arr <- covariate_arr[,,paste0("PC", 1:19)]
all(dimnames(covariate_arr)[[3]] == paste0("PC", 1:19))

#-------------------------------------------------------------------------------
#### Write out object ####
#-------------------------------------------------------------------------------

br_shp <- br %>% rename(codigo_ibg = code_muni)

save(
  chagas_arr,
  pop_arr,
  br_shp,
  covariate_arr,
  file = "chagas_data.Rdata"
)

rm(list=ls())
load("chagas_data.Rdata")


# The remainder of this script is OBSOLETE.
stop("Script executed correctly but the remainder of this script is OBSOLETE")

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









