# Chagas Analyis
setwd("~/chagas_modeling")
library(tidyverse)
library(sf)
library(reshape2)
library(tmap)
tmap_mode("view")
library(posterior)
library(cmdstanr)
library(geobr)
library(cowplot)



load("chagas_data.Rdata")
load("./from_ayesha/cleaned_data/pop_all.Rdata")

to_drop <- br_shp %>% slice(1524, 3498, 5564) %>% pull(muni_code)
br_shp <- br_shp %>%filter(!muni_code %in% to_drop)

# Get map of UFs
ufs <- geobr::read_state()

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(x))
}
# chagas_sample <- readRDS("mcmc_out/chagas_offset.RDS")

# Quick analysis: goodness-of-fit 
chagas_count_v <- as.vector(chagas_offset_count)
pois_ll <- function(p) {
  sum(-dpois(chagas_count_v, p, log=TRUE))
}
pois_ll <- optim(c(0.01), pois_ll, method = "L-BFGS-B", lower = 2e-16, upper = 10)

zip_ll <- function(x) {
  pi <- x[1]
  lambda <- x[2]
  ll = 0
  
  for (i in chagas_count_v){
    if (i == 0){
      ll = ll + -log(pi + (1-pi)*dpois(0,lambda))
    }
    else {
      ll = ll + -log(1-pi) + -(dpois(i, lambda, log=TRUE))
    }
  }
  return(ll)
}
zip_ll <- optim(c(0.9,0.01), zip_ll, method = "L-BFGS-B", lower = c(0, 0.0000001), upper = c(0.99999, 10))

lr <- -2*(zip_ll$value - pois_ll$value)
dchisq(lr, length(chagas_count_v) - 1)
v <- 1:200000
ggplot() + 
  geom_line(aes(v, dchisq(v,length(chagas_count_v) - 1 ))) + 
  geom_vline(aes(xintercept =  lr))


# Analysis of time trend
chagas_offset_count %>% group_by(year) %>% summarize(count= sum(count)) %>%
  ggplot() + geom_point(aes(x=year, y=count))


load <- TRUE

# Predicted counts
if (load) {
  load("E_y_summary.Rdata")
} else {
  E_y1 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-1-4e0ded.csv",
                             variables = "E_y", format = "draws_list")
  E_y2 <- read_cmdstan_csv(
      "mcmc_out/05_knorr_held-202304192205-2-4e0ded.csv",
                             variables = "E_y", format = "draws_list")
  E_y3 <- read_cmdstan_csv(
      "mcmc_out/05_knorr_held-202304192205-3-4e0ded.csv",
                             variables = "E_y", format = "draws_list")
  E_y4 <- read_cmdstan_csv(
      "mcmc_out/05_knorr_held-202304192205-4-4e0ded.csv",
                             variables = "E_y", format = "draws_list")
  
  E_y <- bind_draws(x = E_y1$post_warmup_draws, E_y2$post_warmup_draws, 
                    E_y3$post_warmup_draws, E_y4$post_warmup_draws, along = "chain")
  rm(E_y1, E_y2, E_y3, E_y4)
  gc()
  
  # Since the datset is so big need to build the summary table by hand
  E_y_names <- variables(E_y)
  
  start.time <- Sys.time()
  E_y_summary <- map_dfr(E_y_names, ~{
    d <- extract_variable(E_y, .)
    tibble(variable = ., median = median(d), rhat = rhat(d), ess_median = ess_median(d))
  })
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken ## About 23 minutes
  
  save(E_y_summary, file="E_y_summary.Rdata")
}

# phi, alpha, theta, mu
if (load) {
  load("phi_alpha_summary.Rdata")
} else {
  
  # Repeat for Phi, alpha, and theta, which are much samller 
  phi_alpha1 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-1-4e0ded.csv",
    variables = c("phi_pi","alpha_pi", "theta_pi", "mu_pi"), format = "draws_list")
  phi_alpha2 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-2-4e0ded.csv",
    variables =  c("phi_pi","alpha_pi", "theta_pi", "mu_pi"), format = "draws_list")
  phi_alpha3 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-3-4e0ded.csv",
    variables =  c("phi_pi","alpha_pi", "theta_pi", "mu_pi"), format = "draws_list")
  phi_alpha4 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-4-4e0ded.csv",
    variables =  c("phi_pi","alpha_pi", "theta_pi", "mu_pi"), format = "draws_list")
  
  phi_alpha <- bind_draws(x = phi_alpha1$post_warmup_draws, phi_alpha2$post_warmup_draws, 
                    phi_alpha3$post_warmup_draws, phi_alpha4$post_warmup_draws, along = "chain")
  
  phi_alpha_names <- variables(phi_alpha)
  phi_alpha_summary <- map_dfr(phi_alpha_names, ~{
    d <- extract_variable(phi_alpha, .)
    tibble(variable = ., median = median(d), rhat = rhat(d), ess_median = ess_median(d))
  })
  
  save(phi_alpha_summary, file= "phi_alpha_summary.Rdata")

}

# Standard deviations 
if (load) {
  load("sigma_summary.Rdata")
} else {
  sigma_1 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-1-4e0ded.csv",
    variables = c("sigma_alpha_pi", "sigma_phi_pi", "sigma_theta_pi"), format = "draws_list")
  sigma_2 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-2-4e0ded.csv",
    variables =  c("sigma_alpha_pi", "sigma_phi_pi", "sigma_theta_pi"), format = "draws_list")
  sigma_3 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-3-4e0ded.csv",
    variables = c("sigma_alpha_pi", "sigma_phi_pi", "sigma_theta_pi"), format = "draws_list")
  sigma_4 <- read_cmdstan_csv(
    "mcmc_out/05_knorr_held-202304192205-4-4e0ded.csv",
    variables = c("sigma_alpha_pi", "sigma_phi_pi", "sigma_theta_pi"), format = "draws_list")
  
  sigma <- bind_draws(x = sigma_1$post_warmup_draws, sigma_2$post_warmup_draws, 
                          sigma_3$post_warmup_draws, sigma_4$post_warmup_draws, along = "chain")
  
  sigma_names <- variables(sigma)
  sigma_summary <- map_dfr(sigma_names, ~{
    d <- extract_variable(sigma, .)
    tibble(variable = ., median = median(d), rhat = rhat(d), ess_bulk= ess_bulk(d))
  })
  
  save(sigma_summary, file= "sigma_summary.Rdata")
}

bayesplot::mcmc_trace(sigma)
bayesplot::mcmc_dens(sigma)

mu <- phi_alpha_summary %>%
  filter(str_detect(variable, "mu"))
mu <- mu %>% pull(median)

sigma_phi <- sigma_summary %>% filter(str_detect(variable, "sigma_phi_pi")) %>% pull(median)
sigma_theta <- sigma_summary %>% filter(str_detect(variable, "sigma_theta_pi")) %>% pull(median)
sigma_alpha <- sigma_summary %>% filter(str_detect(variable, "sigma_alpha_pi")) %>% pull(median)

# Overall spatial pattern
overall_shp <- br_shp %>%
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "phi_pi")) %>%
      mutate(phi = median * sigma_phi)  %>%
      select(phi)
  ) %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "theta_pi"))  %>%
      mutate(theta = median * sigma_theta) %>%
      select(theta)
  ) %>%
  mutate(prob = 1-sigmoid(mu + theta+ phi))
  
overall_plot <- tm_shape(overall_shp) + 
  tm_fill(col = "prob", style="cont", title="Pr")  +
  tm_shape(ufs) +
  tm_borders(col="white") + 
  tm_text("abbrev_state", size = 0.5, col = "white") +
  tm_layout(main.title= 'Pr(count>0): Overall Spatial Pattern', 
            main.title.size = 1)

tmap_save(overall_plot, filename = "fig1.png", width=7, height = 5)

# Plot phi and theta
phi_shp <- br_shp %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "phi_pi")) %>%
      mutate(median = median * sigma_phi) %>%
      mutate(median = 1-sigmoid(median))
  )
theta_shp <- br_shp %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "theta_pi"))  %>%
      mutate(median= median * sigma_theta) %>%
      mutate(median =  1-sigmoid(median))
  )

br_phi_theta <- tmap_arrange(
  tm_shape(phi_shp) + 
    tm_fill(col = "median", style="cont", title="Phi")  +
    tm_shape(ufs) +
    tm_borders() + 
    tm_text("abbrev_state", size = 0.5) +
    tm_layout(main.title= 'Pr(count>0): Clustering term', 
              main.title.size = 1) , 
  tm_shape(theta_shp) + 
    tm_fill(col = "median", style="cont", title="Theta") +
    tm_shape(ufs) +
    tm_borders() + 
    tm_text("abbrev_state", size = 0.5) +
    tm_layout(main.title= 'Pr(count>0):  Heterogeneity Term',
              main.title.size = 1),
  nrow=1
)

# tmap_mode("plot")
# br_phi_theta

tmap_save(br_phi_theta, filename = "fig2.png", width=7, height = 5)

# Figure 1
fig1_df <- phi_alpha_summary %>% 
  select(variable, median) %>%
  filter(!str_detect(variable, "alpha|mu")) %>%
  mutate(
    variable = str_replace_all(variable, "\\[", "_"),
    variable = str_replace_all(variable, "\\[|\\]|_pi", "")
    ) %>%
  separate(variable, into = c("variable", "municipality")) %>% 
  pivot_wider(names_from = variable, values_from = median) 

fig1_df <- bind_cols(fig1_df, mu = mu)

fig1_df <- fig1_df %>%
  mutate(pi = theta+phi) %>%
  select(municipality, pi) # %>%
  # mutate(pi = sigmoid(pi))


fig1_df <- bind_cols(br_shp, fig1_df)

tm_shape(fig1_df) +
  tm_fill(col="pi")


# Group E_y by state
E_y_summary <- E_y_summary %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,", ".")) %>%
  separate(variable, sep="\\.", into = c("variable", "year", "muni"))

E_y_shp <- E_y_summary %>%
  select(year, muni, median) %>%
  pivot_wider(names_from = year, values_from = median, names_prefix = "year_")

E_y_shp <- bind_cols(br_shp, E_y_shp)

E_y_shp <- E_y_shp %>%
  rowwise() %>%
  mutate(populacao = as.numeric(populacao)) %>%
  mutate(across(contains("year_"), ~./populacao)) %>%
  mutate(sum_E = sum(c_across(starts_with("year"))))

# Calculate state sum
E_y_shp %>%
  st_drop_geometry() %>%
  group_by(uf) %>% 
  summarize(total_E_y = sum(sum_E)) %>%
  arrange(-total_E_y)

# Possible that theta and phi are not identified, or that theta is dominating
# over phi. Possible that BYM2-convolution can help...
E_y_shp %>%
  st_drop_geometry() %>%
  filter(uf == "MG") %>%
  pivot_longer(contains("year"), names_to = "year") %>%
  mutate(year = as.numeric(str_replace(year, "year_", ""))) %>%
  ggplot() + 
    geom_line(aes(x = year, y = value, group=nome))


# Alpha over time 
alpha_plot <- phi_alpha_summary %>%
  filter(str_detect(variable, "alpha")) %>%
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,", ".")) %>%
  separate(variable, sep="\\.", into = c("variable", "year")) %>%
  mutate(year = as.numeric(year)) %>%
  ggplot() + 
    geom_line(aes(2000:2018, sigmoid(mu + sigma_alpha*median))) + 
    xlab("Year")  + 
    ylab("Probability") + 
    theme_cowplot()

alpha_plot

# Relationships with naive IR









ggplot() + geom_point(aes(
  phi_shp$ess_median,
  theta_shp$ess_median
))

# Check pairs plot
ggplot() + geom_point(
  aes(
    extract_variable(phi_alpha, "phi_pi[1000]"),
    extract_variable(phi_alpha, "theta_pi[1000]")
  )
)

# Check convergence
ggplot() + geom_histogram(aes(E_y_summary$rhat)) + geom_vline(aes(xintercept = 1.01))
ggplot() + geom_histogram(aes(phi_alpha_summary$rhat)) + geom_vline(aes(xintercept = 1.01))

mean(E_y$ess_median)
mean(phi_alpha_summary$ess_median)

# Some variables clearly didn't converge:
E_y_summary %>%
  arrange(-rhat) %>% 
  slice_head(prop=0.01)

E_y_summary %>%
  arrange(ess_bulk) %>% 
  slice_head(prop=0.01)

chagas_summary %>% summarize(rhat_below_1.01 = mean(rhat < 1.01),
                             rhat_below_1.1 = mean(rhat < 1.1),
                             mean_ESS = mean(ess_median)) 


# Pick a few municipalities and plot their parameters over time
chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median))

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,190\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) 

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median))

# Test mapping
IR <- chagas_summary %>%
  filter(str_detect(variable, "^E_y\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")
IR <- IR %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

# count <- as.data.frame(t(chagas_offset_count[1:5, ])) 
# colnames(count) = c("MLE_y1","MLE_y2", "MLE_y3", "MLE_y4", "MLE_y5")
# 
# IR <- bind_cols(IR, count) 
IR <- bind_cols(br_shp, IR)

IR <- IR %>% 
  mutate(populacao = log(as.numeric(populacao))) %>%
  mutate(across(contains("year"), ~./populacao))


tm_shape(IR) + tm_polygons(col = "year4", style = "cont", id = "nome")





pi <- chagas_summary %>%
  filter(str_detect(variable, "^pi\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")

pi <- pi %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

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