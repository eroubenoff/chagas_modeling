# Chagas Analyis
setwd("~/chagas_modeling")
library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)
library(tmap)
tmap_mode("view")

chagas_sample <- readRDS("mcmc_out/chagas_offset.RDS")

## Summarize the parameters separately 
E_y <- chagas_sample$draws("E_y", format = "draws_list") %>% posterior::summarize_draws(median, rhat, ess_bulk)

chagas_summary <- chagas_sample$summary()
chagas_summary
# Check convergence
ggplot() + geom_histogram(aes(chagas_summary$rhat)) + geom_vline(aes(xintercept = 1.01))

mean(chagas_summary$ess_bulk)
sd(chagas_summary$ess_bulk)

# Some variables clearly didn't converge:
chagas_summary %>%
  arrange(-rhat) %>% 
  slice_head(prop=0.01)

chagas_summary %>%
  arrange(ess_bulk) %>% 
  slice_head(prop=0.01)

# I think it is mostly phi_lambda when counts are low
chagas_summary %>% filter(str_detect(variable, "psi"))


# Look at some trace plots
mcmc_areas(chagas_sample$draws("beta0"))
mcmc_areas(chagas_sample$draws("alpha0"))
mcmc_areas(chagas_sample$draws("rho_pi"))
mcmc_areas(chagas_sample$draws("sigma_pi"))
# mcmc_areas(chagas_sample$draws("rho_lambda"))


# Pick a few municipalities and plot their parameters over time
chagas_summary %>% 
  filter(str_detect(variable, "^pi\\[[0-9]*,12\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) + geom_errorbar(aes(x = t,ymin = q5, ymax = q95))

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,190\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) + geom_errorbar(aes(x = t,ymin = q5, ymax = q95))


chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) + geom_errorbar(aes(x = t,ymin = q5, ymax = q95))

# Test mapping
IR <- chagas_summary %>%
  filter(str_detect(variable, "^E_y\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")
IR <- IR %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

count <- as.data.frame(t(chagas_offset_count[1:5, ])) 
colnames(count) = c("MLE_y1","MLE_y2", "MLE_y3", "MLE_y4", "MLE_y5")

IR <- bind_cols(IR, count) 
IR <- bind_cols(br_shp, IR)

IR <- IR %>% 
  mutate(populacao = as.numeric(populacao)) #%>%
# mutate(across(c(year1, year2, year3, year4, year4, MLE_y1, MLE_y2, MLE_y3, MLE_y4, MLE_y5), ~./populacao))


tmap_arrange(
  tm_shape(IR) + tm_polygons(col = "year4", style = "cont"),
  tm_shape(IR) + tm_polygons(col = "MLE_y4", style = "cont")
)

# Plot of IR and MLE
ggplot(IR) +
  geom_point(aes(x = year1, y = MLE_y1)) + 
  geom_point(aes(x = year2, y = MLE_y2)) + 
  geom_point(aes(x = year3, y = MLE_y3)) + 
  geom_point(aes(x = year4, y = MLE_y4)) + 
  geom_point(aes(x = year5, y = MLE_y5)) 





pi <- chagas_summary %>%
  filter(str_detect(variable, "^pi\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")

pi <- pi %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(x))
}
pi <- pi %>%
  mutate(across(starts_with("year"), sigmoid))
pi <- bind_cols(br_shp, pi)
pi <- bind_cols(pi, count)
pi <- pi %>%
  mutate(populacao = as.numeric(populacao)) %>%
  mutate(MLE_y1 = MLE_y1/populacao,
         MLE_y2 = MLE_y2/populacao,
         MLE_y3 = MLE_y3/populacao,
         MLE_y4 = MLE_y4/populacao,
         MLE_y5 = MLE_y5/populacao
  )

tmap_arrange(
  
  tm_shape(pi) + tm_polygons(col = "year4", style = "cont"),
  tm_shape(pi) + tm_polygons(col = "MLE_y4", style = "cont")
)