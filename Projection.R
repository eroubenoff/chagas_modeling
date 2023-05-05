# Projection code
library(tidyverse)
pop_2019 <- readxl::read_excel("Projection/estimativa_dou_2019.xls", 2, skip=1)
pop_2020 <- readxl::read_excel("Projection/estimativa_dou_2020.xls", 2, skip=1)
pop_2021 <- readxl::read_excel("Projection/estimativa_dou_2021.xls", 2, skip=1)

colnames(pop_2019) <- c("UF", "UF_co", "muni_code", "muni_name", "pop_2019")
colnames(pop_2020) <- c("UF", "UF_co", "muni_code", "muni_name", "pop_2020")
colnames(pop_2021) <- c("UF", "UF_co", "muni_code", "muni_name", "pop_2021")

pop_2019 <- pop_2019 %>%
  mutate(pop_2019 = as.numeric(pop_2019)) %>% 
  drop_na() 

pop_2020 <- pop_2020 %>%
  mutate(pop_2020 = as.numeric(pop_2020)) %>%
  drop_na() 

pop_2021 <- pop_2021 %>%
  mutate(pop_2021 = as.numeric(pop_2021)) %>%
  drop_na() 

pop <- full_join(pop_2019, pop_2020) %>% full_join(pop_2021)

# 60 municipalities were not present in all years:
pop %>% filter(is.na(pop_2019) | is.na(pop_2020) | is.na(pop_2021)) %>% View()

# For those with 2019 but not 2020, take the 2019 value:
pop[is.na(pop$pop_2019),]

