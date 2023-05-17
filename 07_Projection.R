#-------------------------------------------------------------------------------
# Spatial offset model (non-covariate), with time, projecting forward
# Spatial model only with offset
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html
# https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/bym2_offset_only.stan
# download.file("https://raw.githubusercontent.com/stan-dev/example-models/master/knitr/car-iar-poisson/nb_data_funs.R",
#               destfile= "nb_data_funs.R")
#-------------------------------------------------------------------------------

setwd("~/chagas_modeling")
library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)
library(tmap)
library(posterior)
tmap_mode("view")

load("chagas_data.Rdata")

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(-x))
}


# Load the mcmc data
main_model <- bind_rows(
  readRDS("phi_alpha_summary.RDS"),
  readRDS("delta_lambda_summary.RDS"),
  readRDS("delta_pi_summary.RDS"))
  
covariate_model <- bind_rows(
  readRDS("covariate_phi_alpha_summary.RDS"),
  readRDS("covariate_delta_lambda_summary.RDS"),
  readRDS("covariate_delta_pi_summary.RDS"))



## Need to load in all of the indicators for population projections
flist <- list.files("br_2010_census_por_municipios", full.names = TRUE)
flist <- flist %>% str_subset("Tabela 1\\.")

br_2010 <- map_dfr(flist, ~{
  readxl::read_excel(., skip = 6) %>% 
    select(name = 1, pop = 2, code_muni = 11)
})

br_2010 <- br_2010 %>% filter(name != "Total")



project_arr <- function(pop_arr, n) {
  
  # Look at the average annual growth rate over the course of history
  n_already <- dim(pop_arr)[1]
  growth_rates <- exp((log(pop_arr[n_already, ]) - log(pop_arr[1,])) / n_already) - 1
  projected_array <- array(0, dim = list(n, dim(pop_arr)[2]))
  dimnames(projected_array) <- list(seq(max(as.numeric(dimnames(pop_arr)[[1]])) + 1, length.out = n),
                                    dimnames(pop_arr)[[2]])
  
  for (i in 1:n) {
    projected_array[i, ] <- round(pop_arr[n_already,] * (1+growth_rates)^i)
  }
  
  return(list(growth_rates, projected_array))
  
}

n_projection_years <- 10
projected <- project_arr(pop_arr, n_projection_years)
growth_rates <- projected[[1]]
projected <- projected[[2]]
pop_projected <- abind::abind(pop_arr, projected, along=1)

br_proj <- project_arr(matrix(rowSums(pop_arr), 
                   ncol = 1, 
                   dimnames = list(dimnames(pop_arr)[[1]], "BR")), 
            n = n_projection_years)


# Plot of overall BR population
ggplot() + 
  geom_line(aes(x = seq(2001, 2029), y = c(rowSums(pop_arr), rowSums(projected)))) +
  geom_vline(aes(xintercept = 2019), color="red") + 
  cowplot::theme_cowplot()

# Histogram of growth rates
growth_rates_hist <- ggplot() + geom_histogram(aes(100*growth_rates)) +
  geom_vline(aes(xintercept = 100*br_proj[[1]]), color = "red") +
  xlab("Rate (%) ") + 
  ylab("") + 
  ggtitle("Municipio Growth Rates, 2001-2019") + 
  cowplot::theme_cowplot()  

ggsave("growth_rates_hist.png", growth_rates_hist, height = 5, width = 5)

# Map of growth rates
tmap_mode("plot")
ufs <- geobr::read_state()
growth_rates_map <-  br_shp %>% bind_cols(growth_rates = 100*growth_rates) %>%
  tm_shape() + 
  tm_fill(col = "growth_rates", title = "Growth Rate, %", style = "cont")  + 
  tm_shape(ufs) +
  tm_borders(col="black") + 
  tm_text("abbrev_state", size = 0.5, col = "black") +
  tm_layout(main.title= 'Annual Growth Rate: 2001-2019, %', 
            main.title.size = 1)

tmap_save(growth_rates_map, "growth_rates_map.png", height = 5, width = 5)

# Projected population change, 2019-2029:
pop_change_map <-  br_shp %>% bind_cols(pop_change = projected[10,] - pop_arr[19,]) %>%
  tm_shape() + 
  tm_fill(col = "pop_change", title = "Increase", style = "cont")+ 
  tm_shape(ufs) +
  tm_borders(col="black") + 
  tm_text("abbrev_state", size = 0.5, col = "black") +
  tm_layout(main.title= 'Projected Population Increase, 2019-2029', 
            main.title.size = 1)

tmap_save(pop_change_map, "pop_change_map.png", height = 5, width = 5)



source("03_chagas_helper.R")
scaling_factor<- create_stan_data(chagas_arr = chagas_arr, pop_arr = pop_arr, br_shp = br_shp, testing = FALSE , covariate_arr=covariate_arr, nvar= 6)$scaling_factor


# Project using the covariate data
project_chagas <- function(n_reps, years, model, projected_pop, scaling_factor, projected_covariates = NULL) {
  mu_lambda <- model %>% filter(str_detect(variable, "^mu_lambda")) %>% pull(median) %>% as.numeric()
  mu_pi <- model %>% filter(str_detect(variable, "^mu_pi")) %>% pull(median) %>% as.numeric()
  
  alpha_lambda <- model %>% filter(str_detect(variable, "^alpha_lambda")) %>% pull(median)  %>% as.numeric()
  sigma_alpha_lambda <- model %>% filter(str_detect(variable, "^sigma_alpha_lambda")) %>% pull(median)  %>% as.numeric()
  alpha_pi <- model %>% filter(str_detect(variable, "^alpha_pi")) %>% pull(median)  %>% as.numeric()
  sigma_alpha_pi <- model %>% filter(str_detect(variable, "^sigma_alpha_pi")) %>% pull(median)  %>% as.numeric()
  
  if (!is.null(projected_covariates)) {
    beta_lambda <- model %>% filter(str_detect(variable, "^beta_lambda")) %>% pull(median) %>% as.numeric()
    beta_pi <- model %>% filter(str_detect(variable, "^beta_pi")) %>% pull(median)  %>% as.numeric()
  } else {
    beta_lambda <- 0
    beta_pi <- 0
  }
  
  
  theta_lambda <- model %>% filter(str_detect(variable, "^theta_lambda")) %>% pull(median) %>% as.numeric()
  theta_pi <- model %>% filter(str_detect(variable, "^theta_pi")) %>% pull(median) %>% as.numeric()
  
  phi_pi <- model %>% filter(str_detect(variable, "^phi_pi")) %>% pull(median) %>% as.numeric()
  rho_pi <- model %>% filter(str_detect(variable, "^rho_pi")) %>% pull(median) %>% as.numeric()
  
  sigma_convolved_lambda <- model %>% filter(str_detect(variable, "^sigma_convolved_lambda")) %>% pull(median)  %>% as.numeric()
  sigma_convolved_pi <- model %>% filter(str_detect(variable, "^sigma_convolved_pi")) %>% pull(median)  %>% as.numeric()
  
  
  max_years <- 19+years
  
  delta_lambda <- model %>% 
    filter(str_detect(variable, "^delta_lambda")) %>% 
    mutate(year = str_extract(variable, "(?<=\\[)[0-9]*"),
           muni_no = str_extract(variable, "[0-9]*(?=\\])")) %>% 
    mutate(across(c(year, muni_no, median, rhat, ess), ~as.numeric(.))) %>% 
    acast(year~muni_no, value.var = "median") %>%
    abind::abind(array(0, c(years, length(theta_pi))), along =1)
  
  delta_pi <- model %>% 
    filter(str_detect(variable, "^delta_pi")) %>% 
    mutate(year = str_extract(variable, "(?<=\\[)[0-9]*"),
           muni_no = str_extract(variable, "[0-9]*(?=\\])")) %>% 
    mutate(across(c(year, muni_no, median, rhat, ess), ~as.numeric(.))) %>%
    acast(year~muni_no, value.var = "median")  %>%
    abind::abind(array(0, c(years, length(theta_pi))), along =1)
  

  
  
  sigma_delta_lambda <- model %>% filter(str_detect(variable, "^sigma_delta_lambda")) %>% pull(median)  %>% as.numeric()
  sigma_delta_pi <- model %>% filter(str_detect(variable, "^sigma_delta_pi")) %>% pull(median)  %>% as.numeric()
  
  
  ## Assemble model
  pi_f <- function(mu, alpha, sigma_alpha, beta, theta, phi, sigma, delta, sigma_delta, rho, covariates, scaling_factor) {
    
    if (!is.null(covariates)) {
      covs <- as.vector(as.matrix(covariates) %*% beta)
    } else {
      covs <- 0 
    }
    return(sigmoid(mu + alpha + covs + (sqrt(rho/scaling_factor)*phi + sqrt(1 - rho)*theta)*sigma + delta*sigma_delta))
  }
  
  lambda_f <- function(mu, alpha, sigma_alpha, beta, theta, delta, sigma, sigma_delta, covariates) {
    
    if (!is.null(covariates)) {
      covs <- as.vector(as.matrix(covariates) %*% beta)
    } else {
      covs <- 0 
    }
    return(exp(mu + alpha + covs + theta*sigma + delta*sigma_delta))
  }
  
  
  
  ret = vector("list", length = n_reps)
  alpha_lambda_old <- alpha_lambda
  alpha_pi_old <- alpha_pi
  for (n in 1:n_reps) {
    ## Assemble return object
    projected_chagas <- projected_pop[,]
    projected_chagas[,] <- 0
    pi <- projected_chagas
    lambda <- projected_chagas
    E_y <- projected_chagas
  
    ## Simulate random walk in alpha
    alpha_lambda <- alpha_lambda_old
    alpha_pi <- alpha_pi_old
    for (i in 1:years) {
      alpha_lambda <- c(alpha_lambda, rnorm(1, mean = tail(alpha_lambda, 1), sd = sigma_alpha_lambda))
      alpha_pi <- c(alpha_pi, rnorm(1, mean = tail(alpha_pi, 1), sd = sigma_alpha_pi))
    }
    
    
    # for (i in 1:ncol(delta_lambda)) {
      # m1 <- delta_lambda[1:19, i]
      # delta_lambda[20:max_years, i] <- sample(m1, size=years, replace = TRUE)
      # m1 <- delta_pi[1:19, i]
      # delta_pi[20:max_years, i] <- sample(m1, size=years, replace = TRUE)
    # }
    
    ## Do projected 
    for (i in 1:max_years) {
      pi[i,] <- pi_f(mu = mu_pi, alpha = alpha_pi[i], sigma_alpha = sigma_alpha_pi,
                     beta = beta_pi, theta = theta_pi, phi = phi_pi, sigma = sigma_convolved_pi,
                     delta = delta_pi[i,], 
                     sigma_delta = sigma_delta_pi, rho = rho_pi, covariates = projected_covariates[i,,],
                     scaling_factor = scaling_factor)
      lambda[i, ] <- lambda_f(mu = mu_lambda, alpha = alpha_lambda[i], sigma_alpha = sigma_alpha_lambda,
                              beta = beta_lambda, theta = theta_lambda, sigma = sigma_convolved_lambda,
                              delta = delta_lambda[i,], 
                              sigma_delta = sigma_delta_lambda, covariates = projected_covariates[i,,])
      
      E_y[i,] <- (1-pi[i,]) * lambda[i,] * projected_pop[i,]
    
    }
    
    ret[[n]] <- E_y
  }
  
  return(abind::abind(ret, along = 3))
  
}


set.seed(494949)

projected_chagas_main <- project_chagas(n_reps = 1000,years = 10, model = main_model, 
                projected_pop = pop_projected, 
                scaling_factor = scaling_factor)


projected_chagas_covariate <-  project_chagas(n_reps = 1000, years = 10, model = covariate_model, 
                  projected_pop = pop_projected, 
                  projected_covariates = covariate_arr[as.character(2001:2030),,1:6],
                  scaling_factor = scaling_factor)

projected_chagas_main <- bind_rows(
  melt(projected_chagas_main[1:19, ,1]), 
  melt(projected_chagas_main[20:29,,])
) %>% tibble() %>% 
  rename(year = 1, muni = 2, value = 3)

projected_chagas_covariate <- bind_rows(
  melt(projected_chagas_covariate[1:19, ,1]), 
  melt(projected_chagas_covariate[20:29,,])
) %>% tibble() %>%
  rename(year = 1, muni = 2, value = 3)




projected_chagas_main_summary <- projected_chagas_main %>%
  group_by(year, Var3) %>% summarize(value = sum(value)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(med= median(value), qlow = quantile(value, 0.25), qupp = quantile(value, 0.75))


projected_chagas_covariate_summary <- projected_chagas_covariate %>% 
  group_by(year, Var3) %>% summarize(value = sum(value)) %>% 
  ungroup() %>% 
  group_by(year) %>% 
  summarize(med = median(value), qlow = quantile(value, 0.25), qupp = quantile(value, 0.75))


sum(tail(rowSums(chagas_arr), 10))

projected_chagas_main_summary %>% 
  filter(year > 2019) %>% 
  summarize(across(c(med, qlow, qupp), sum))

projected_chagas_covariate_summary %>% 
  filter(year > 2019) %>% 
  summarize(across(c(med, qlow, qupp), sum))

br_pop_total <- rowSums(projected)
br_pop_total[10] - br_pop_total[1]


proj_plot <- 
cowplot::plot_grid(
  ggplot() + 
  geom_ribbon(aes(year, ymin = qlow, ymax = qupp), data = projected_chagas_main_summary %>% filter(year > 2019), alpha = 0.5) +
  geom_line(aes(year, med), data = projected_chagas_main_summary %>% filter(year > 2019)) +
  geom_line(aes(x=2001:2019, y = rowSums(chagas_arr)), linetype= "dashed") +
  geom_ribbon(aes(x=2001:2019, ymin = rowSums(chagas_arr),ymax= rowSums(chagas_arr))) +
  cowplot::theme_cowplot() + 
  ylab("Count") + 
  xlab("Year") +
  ggtitle("Projected Annual Counts of Chagas Disease", subtitle = "Main Model") +
  expand_limits(y=0),
  
  ggplot() + 
  geom_ribbon(aes(year, ymin = qlow, ymax = qupp), alpha = 0.5, data = projected_chagas_covariate_summary %>% filter(year > 2019)) +
  geom_line(aes(year, med), data = projected_chagas_covariate_summary %>% filter(year > 2019)) +
  geom_line(aes(x=2001:2019, y = rowSums(chagas_arr)), linetype= "dashed") +
  geom_ribbon(aes(x=2001:2019, ymin = rowSums(chagas_arr),ymax= rowSums(chagas_arr))) +
  cowplot::theme_cowplot() + 
  ylab("Count") + 
  xlab("Year") +
  ggtitle("", subtitle = "Climate Covariate Model") +
  expand_limits(y=0),
  ncol =1
  )


ggsave("proj_plot.png", proj_plot, width=7, height=7)

br_shp <- br_shp %>% 
  rowwise() %>%
  mutate(chagas_total = sum(across(paste0("chagas_", 2010:2019)))) %>% 
  left_join(
    projected_chagas_main %>% 
      filter(year > 2019) %>% 
      group_by(muni_code = muni, year) %>% 
      summarize(chagas_proj = median(value)) %>% 
      ungroup(year) %>% 
      summarize(chagas_proj_main = sum(chagas_proj))
  ) %>% 
  left_join(
  projected_chagas_covariate %>% 
    filter(year > 2019) %>% 
    group_by(muni_code = muni, year) %>% 
    summarize(chagas_proj = median(value)) %>% 
    ungroup(year) %>% 
    summarize(chagas_proj_covariate = sum(chagas_proj))
  )

pct_increase <- function(x, y) {
  if (x == 0 & y == 0) {return(NA)}
  else if (x == 0) {return(0)}
  else if (y == 0) {return(0)}
  else {return((x - y)/y*100)}
}

br_shp <- br_shp %>%
  mutate(main_increase = pct_increase(chagas_proj_main, chagas_total),
         covariate_increase = pct_increase(chagas_proj_covariate, chagas_total))

proj_map <- tmap_arrange(
  tm_shape(br_shp) + 
    tm_fill(col = "chagas_proj_main", style = "cont", title = "Count") + 
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'Main Model: Projected Cases 2020-2030', main.title.size = 1),
  tm_shape(br_shp) + 
    tm_fill(col = "main_increase", style = "cont", title = "%")+ 
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black")+
    tm_layout(main.title= 'Main Model: % Increase 2020-2030', main.title.size = 1),
  tm_shape(br_shp) + 
    tm_fill(col = "chagas_proj_covariate", style = "cont", title = "Count") + 
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black")+
    tm_layout(main.title= 'Covariate Model: Projected Cases 2020-2030', main.title.size = 1),
  tm_shape(br_shp) + 
    tm_fill(col = "covariate_increase", style = "cont", title = "%") + 
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black")+
    tm_layout(main.title= 'Covariate Model: % Increase 2020-2030', main.title.size = 1)
)
  
tmap_save(proj_map, "proj_map.png", width = 8, height = 8)

proj_muni <- br_shp %>% 
  st_drop_geometry() %>% 
  ungroup() %>%
  arrange(-chagas_proj_main) %>% 
  select(name_muni, abbrev_state, pop_2019, chagas_total, contains("proj"), contains("increase")) %>% 
  # mutate(chagas_total = percent_rank(chagas_total),
  #        chagas_proj_main = percent_rank(chagas_proj_main),
  #        chagas_proj_covariate= percent_rank(chagas_proj_covariate)
  #        ) %>% 
  filter(chagas_total > 0 , chagas_proj_main > 1 , chagas_proj_covariate > 1) 

proj_muni_main <- proj_muni %>%
    ggplot() + 
      geom_text(aes(chagas_total, chagas_proj_main, label = name_muni)) + 
      xlab("Total Chagas Counts, 2001-2019") + 
      ylab("Predicted Chagas Counts, 2020-2030") + 
    scale_x_log10() + 
    scale_y_log10() +
    cowplot::theme_cowplot()
proj_muni_cov <-  proj_muni %>%
    ggplot() + 
      geom_text(aes(chagas_total, chagas_proj_covariate, label = name_muni))  +
      xlab("Total Chagas Counts, 2001-2019") + 
      ylab("Predicted Chagas Counts, 2020-2030") + 
    scale_x_log10() + 
    scale_y_log10() +
    cowplot::theme_cowplot()
  
ggsave("proj_muni_main.png", proj_muni_main, height = 7, width = 7)
ggsave("proj_muni_cov.png", proj_muni_cov, height = 7, width = 7)

br_shp %>% 
  st_drop_geometry() %>% 
  arrange(-chagas_proj_main) %>% 
  select(name_muni, abbrev_state, pop_2019, contains("proj"), contains("increase"))


br_shp %>% 
  st_drop_geometry() %>% 
  arrange(-main_increase) %>% 
  select(name_muni, abbrev_state, pop_2019, contains("proj"), contains("increase"))

br_shp %>% 
  st_drop_geometry() %>% 
  arrange(-covariate_increase) %>% 
  select(name_muni, abbrev_state, pop_2019, contains("proj"), contains("increase"))



