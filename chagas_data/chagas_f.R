# Produce the data needed from the tabnet microdata

library(tidyverse)
library(read.dbc)
library(tabulizer)
library(geobr)
# library(poly)
setwd("/nobackup/90days/eroubenoff/chagas/chagas_data")

flist <- list.files(path = "./microdata_download", full.names = TRUE)

chagas <- map_dfr(.x = flist, .f = ~read.dbc(., as.is = TRUE))
str(chagas)

chagas <- chagas %>%
  select(notification_date = DT_NOTIFIC, # dd/mm/aaaa
         uf_of_notification = SG_UF_NOT,
         municipality_of_notification = ID_MUNICIP,
         unidade_of_notification = ID_UNIDADE, # unclear to me how this field differs from SF_UF_NOT
         date_of_first_symptoms = DT_SIN_PRI,
         date_of_birth= DT_NASC,
         age = NU_IDADE_N,  # this column is confusing? 
         sex = CS_SEXO,
         race = CS_RACA, #1: branca 2: preta 3: amarela 4: parda 5: indigena 9: ignorado
         education = CS_ESCOL_N,
         uf_of_residence = SG_UF,
         municipality_of_residence = ID_MN_RESI,
         country = ID_PAIS, # If not resident of brazil
         occupation = ID_OCUPA_N,
         
         triatomines_in_home = PRESENCA, # 1: yes 2: no 3: unknown 9: ignored
         date_of_parasite_encounter = PARASITO,
         use_of_blood_products = HISTORIA, # did they use blood products in last 90 days
         works_with_t_cruzi = MANIPULA, # does the patient work with materiales that have t crusi
         mother_infected = MAECHAGA, # if patient under 9 months, was the mother infected
         uf_of_infection = COUFINF, 
         municipality_of_infection = COMUNINF 
         # There are some other columns but I think this is all we need for now
         )


# get population data and shapefile (from ayesha!)
load("../from_ayesha/cleaned_data/pop_all.Rdata")
pop <- pop.all
br <- st_read("../from_ayesha/municipios_2010.shp", "municipios_2010")


chagas_monthly <- chagas %>%
  mutate(year_of_first_symptoms = lubridate::year(date_of_first_symptoms),
         month_of_first_symptoms = lubridate::month(date_of_first_symptoms)) %>%
  group_by(municipality_of_residence, year_of_first_symptoms, month_of_first_symptoms) %>% 
  summarize(n = n())

chagas_summary <- chagas_monthly %>%
  group_by(municipality_of_residence) %>%
  summarize(Total = sum(n))



br <- br %>% mutate(code_muni_trunc = as.numeric(substr(codigo_ibg, 1, 6)))
summary(as.numeric(chagas_summary$municipality_of_residence) %in% as.numeric(br$code_muni_trunc))
chagas_summary <- chagas_summary %>% mutate_at(vars(municipality_of_residence), ~as.numeric(.))

chagas_shp <- left_join(br, chagas_summary, by = c("code_muni_trunc" = "municipality_of_residence"))

pop_2015 <- pop %>% 
              filter(year == 2015) %>%
              select(code, pop = ps)
            
chagas_shp <- left_join(chagas_shp, pop_2015, by = c("code_muni_trunc"= "code"))

chagas_shp <- chagas_shp %>%
  mutate(Total = if_else(is.na(Total), 0, as.numeric(Total)),
         FOI = Total/as.numeric(pop))


tm_shape(chagas_shp) +
  tm_fill(col = "FOI", style = "jenks")



# Seasonality
chagas_monthly <- chagas_monthly %>%
  group_by(month_of_first_symptoms, year_of_first_symptoms) %>% 
  summarize(n = sum(n))
  

ggplot(chagas_monthly) +
  geom_line(aes(x = month_of_first_symptoms, y = n, col = as.factor(year_of_first_symptoms)))








flist <- list.files(path = "./population_tcu", full.names = TRUE)
for (f in flist) {
  unzip(f, exdir = "./population_tcu/")
}
flist <- list.files("./population_tcu", full.names = TRUE) %>% str_subset("(.dbf|.DBF)")

# In one data source there is 6-digit municipality codes and in the other there are 7-digit
# codes. This thread explains it a little. 
# https://forum.ipums.org/t/municipalities-in-1991-brazilian-census/2847


chagas_summary <- chagas_summary %>%
  pivot_wider(names_from = year_of_first_symptoms, values_from = n, values_fill = 0,
              names_prefix = "n_")

chagas_summary[, "cumulative"] <- rowSums(chagas_summary[, 2:14])

chagas_summary <- chagas_summary %>%
  mutate(municipality_of_residence = as.numeric(municipality_of_residence))

nDigits <- function(x) nchar( trunc( abs(x) ) )
nDigits <- Vectorize(nDigits)
nDigits(chagas_summary$municipality_of_residence) %>% summary
nDigits(br$code_muni) %>% summary

br <- br %>% mutate(code_muni_trunc = as.numeric(substr(code_muni, 1, 6)))

as.numeric(chagas_summary$municipality_of_residence) %in% as.numeric(br$code_muni_trunc)

chagas_shp <- left_join(br, chagas_summary, by = c("code_muni_trunc" = "municipality_of_residence"))

# chagas_shp[chagas_shp$cumulative == 5054, "cumulative"] <- 0
tm_shape(chagas_shp) +
  tm_fill(col = "cumulative", style = "log10")


# Crude FOI
pop <- pop %>% 
  filter(ANO ==2018) %>% 
  select(-ANO) %>% 
  mutate(MUNIC_RES = as.numeric(substr(MUNIC_RES, 1, 6)))



chagas_shp <- chagas_shp %>%
  left_join(pop, by = c("code_muni_trunc" = "MUNIC_RES"))

chagas_shp <- chagas_shp %>%
  mutate(FOI = cumulative/POPULACAO,
         FOI_log = -log10(FOI))

tm_shape(chagas_shp) + 
  tm_fill(col = "FOI", style = "cont")



# Something about these maps isn't looking right... 
l <- list.files("from_tabnet", full.names = TRUE)
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

save(chagas_municipio, file = "chagas_municipio.Rdata")

chagas_municipio <- chagas_municipio %>%
  select(-year, -level, -`UF de residência`, -`Ign/Em Branco` ) %>%
  group_by(`Município de residência`) %>%
  summarize(Total = sum(Total))
  

chagas_municipio <- chagas_municipio %>%
  mutate(muni_code = as.numeric(str_sub(`Município de residência`, 1, 6)))
  
 
chagas_municipio$muni_code %in% br$code_muni_trunc

br_2 <- st_read("../from_ayesha/municipios_2010.shp", "municipios_2010")
br_2 <- br_2 %>%
  mutate(muni_code = as.numeric(substr(codigo_ibg, 1, 6)))

chagas_shp <- left_join(br_2, chagas_municipio, by = c("muni_code" = "muni_code"))

chagas_shp$Total <- replace_na(chagas_shp$Total, 0)
chagas_shp <- chagas_shp %>% mutate(Total = replace_na(Total, 0),
                     FOI = Total / as.numeric(populacao))

chagas_shp[which.max(chagas_shp$FOI),]

tm_shape(chagas_shp) +
  tm_fill(col = "FOI")

























