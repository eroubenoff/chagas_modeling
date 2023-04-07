# Function to read in the Chagas data from tabnet (not the indiviudal data)
# and print some descriptives

library(tidyverse)
library(read.dbc)
library(tabulizer)
library(sf)
library(tmap)
rm(list = ls())
setwd("chagas_modeling/chagas_data")

l <- list.files("from_tabnet", full.names = TRUE)
chagas_tabnet <- map_dfr(l, ~{
  df <- read_delim(., delim = ";", na = "-")
  year <- str_sub(., -8, -5)
  level <- if_else(str_detect(., "UF"), "UF", "municipio")
  df <- df %>%
    mutate(year = year,
           level = level) 
})

chagas_municipio <-  chagas_tabnet %>%
  filter(`Município de residência` != "Total") %>%
  filter(level == "municipio")

chagas_municipio <- chagas_municipio %>%
  select(-year, -level, -`UF de residência`, -`Ign/Em Branco` ) %>%
  group_by(`Município de residência`) %>%
  summarize(Total = sum(Total))


chagas_municipio <- chagas_municipio %>%
  mutate(muni_code = as.numeric(str_sub(`Município de residência`, 1, 6)))



br <- st_read("../from_ayesha/municipios_2010.shp", "municipios_2010")
br <- br %>%
  mutate(muni_code = as.numeric(substr(codigo_ibg, 1, 6)))

as.numeric(chagas_municipio$muni_code) %in% as.numeric(br$muni_code)


chagas_shp <- left_join(br, chagas_municipio, by = c("muni_code" = "muni_code"))

chagas_shp$Total <- replace_na(chagas_shp$Total, 0)
chagas_shp <- chagas_shp %>% mutate(FOI = Total / as.numeric(as.character(populacao)),
                                    FOI_log10 = if_else(FOI == 0, 0, -log10(FOI)))

chagas_shp[which.max(chagas_shp$FOI),]

(tm_shape(chagas_shp) +
  tm_fill(col = "FOI", style = "cont", palette = "Reds")) %>%
  tmap_save(filename = "chagas_plot.png", width =4, height = 4, units="in" )



# Seasonality 

chagas_monthly <- chagas_tabnet %>%
  filter(level == "UF") %>%
  select(-`Município de residência`, -Total, -level, -`Ign/Em Branco`, -`UF de residência`) 

colnames(chagas_monthly) <- c(as.character(1:12), "year")

chagas_monthly <- chagas_monthly %>%
  group_by(year) %>% 
  summarize_all(~sum(., na.rm = TRUE))

chagas_monthly <- chagas_monthly %>% pivot_longer(-year, names_to = "month")

ggplot(chagas_monthly) +
  geom_line(aes(x=month, y = value, group = year)) + 
  stat_smooth(aes(x = month, y = value))


chagas_monthly <- chagas_monthly %>%
  select(-year) %>%
  summarize_all(~sum(., na.rm = TRUE))
  
chagas_monthly <- chagas_monthly %>%
  pivot_longer(1:12, names_to = "month", values_to = "n")



  
(ggplot(chagas_monthly) +
  geom_line(aes(x = as.numeric(month), y= n), size = 2 )) %>%
  ggsave(filename = "chagas_monthly.png", width = 4, height=4, units ="in")
  
  
  
  
  
  
  


