#-------------------------------------------------------------------------------
# Spatial offset model (non-covariate), with time
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
# we recommend running this is a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
# install_cmdstan(cores = 4, version = "2.31.0")
set_cmdstan_path("/Users/eroubenoff/.cmdstan/cmdstan-2.31.0")

source("chagas_helper.R")

load("chagas_data.Rdata")

options("cmdstanr_verbose" = TRUE)

testing <- TRUE 
stan_data <- create_stan_data(chagas_arr, pop_arr, br_shp, testing = testing, covariate = FALSE)
chagas_offset <- cmdstan_model("knorr_held_convolved_both_parts.stan",
                               force_recompile = FALSE)

message("Sampling began at", Sys.time())
chagas_sample <- chagas_offset$sample(data = stan_data,
                      chains=ifelse(testing, 2, 4),
                      parallel_chains=4,
                      iter_warmup=ifelse(testing, 500, 2000),
                      iter_sampling=ifelse(testing,100, 1000),
                      output_dir = "mcmc_out",
                      save_latent_dynamics = TRUE,
                      adapt_delta = 0.8,
                      refresh=25,
                      save_warmup=FALSE)


stop()

chagas_sample <- as_cmdstan_fit(c(
  "mcmc_out/knorr_held_convolved_both_parts-diagnostic-202304271505-3-53a251.csv",
  "mcmc_out/knorr_held_convolved_both_parts-diagnostic-202304271505-4-53a251.csv"
))
chagas_summary <- chagas_sample$summary() 
warnings()
chagas_summary

mu <- chagas_summary %>% filter(str_detect(variable, "mu")) 
mu
alpha <- chagas_summary %>% filter(str_detect(variable, "alpha")) 
alpha

phi_pi <- chagas_summary %>% filter(str_detect(variable, "phi_pi")) %>% pull(median)
phi_lambda <- chagas_summary %>% filter(str_detect(variable, "phi_lambda")) %>% pull(median)
ggplot() + 
  geom_histogram(aes(phi_pi, fill = "phi pi")) + 
  geom_histogram(aes(phi_lambda, fill = "phi lambda"))

theta_pi <- chagas_summary %>% filter(str_detect(variable, "theta_pi")) %>% pull(median)
theta_lambda <- chagas_summary %>% filter(str_detect(variable, "theta_lambda")) %>% pull(median)
ggplot() + 
  geom_histogram(aes(theta_pi, fill = "theta pi")) + 
  geom_histogram(aes(theta_lambda, fill = "theta lambda"))

rho <- chagas_summary %>% filter(str_detect(variable, "rho")) 
rho
sigma <- chagas_summary %>% filter(str_detect(variable, "sigma_convolved")) 
sigma

# Total spatial term:
mu_pi <- mu %>% filter(variable == "mu_pi") %>% pull(median) 
mu_lambda <- mu %>% filter(variable == "mu_lambda") %>% pull(median) 
rho_pi <- rho %>% filter(variable == "rho_pi") %>% pull(median)
rho_lambda <- rho %>% filter(variable == "rho_lambda") %>% pull(median)
sigma_pi <- sigma %>% filter(variable == "sigma_convolved_pi") %>% pull(median)
sigma_lambda <- sigma %>% filter(variable == "sigma_convolved_lambda") %>% pull(median)

pi_spatial <- mu_pi + sigma_pi*(sqrt(rho_pi/stan_data$scaling_factor)*phi_pi + sqrt(1-rho_pi)*theta_pi)
pi_spatial <- sigmoid(pi_spatial)
hist(pi_spatial)

lambda_spatial <- mu_lambda + sigma_lambda*(sqrt(rho_lambda/stan_data$scaling_factor)*phi_lambda + sqrt(1-rho_lambda)*theta_lambda)
lambda_spatial <- exp(lambda_spatial)
hist(lambda_spatial)

# Some diagnostics
# Rhat
chagas_summary %>% arrange(-rhat) 
chagas_summary %>%
  filter(str_detect(variable, "rho"))

# % under 1.1 and 1.01
chagas_summary %>% summarize(rhat_below_1.01 = mean(rhat < 1.01),
                             rhat_below_1.1 = mean(rhat < 1.1))
chagas_summary %>% arrange(ess_bulk) 
chagas_summary %>% ggplot() + geom_histogram(aes(rhat))

## Traceplot of all single or vector parameters 
posterior_cp <- as_draws_array(chagas_sample$draws())
np_cp <- bayesplot::nuts_params(chagas_sample)
mcmc_trace(chagas_sample$draws("mu_pi"))
mcmc_trace(chagas_sample$draws("mu_lambda"))

mcmc_trace(chagas_sample$draws(c("alpha_pi", "sigma_alpha_pi")))
mcmc_trace(chagas_sample$draws(c("alpha_lambda", "sigma_alpha_lambda")))

mcmc_trace(chagas_sample$draws(c("rho_pi", "rho_lambda")))
mcmc_pairs(posterior_cp, pars = c("rho_pi", "rho_lambda"), np=np_cp)

mcmc_trace(chagas_sample$draws(c("sigma_convolved_pi", "sigma_convolved_lambda")))
mcmc_pairs(posterior_cp, pars = c("sigma_convolved_pi", "sigma_convolved_lambda"), np=np_cp)

mcmc_trace(chagas_sample$draws(c("sigma_delta_pi", "sigma_delta_lambda")))
mcmc_pairs(posterior_cp, pars = c("sigma_delta_pi", "sigma_delta_lambda"), np=np_cp)

## Now to do the large vectors of params, inclu. theta phi and delta
mcmc_trace(chagas_sample$draws("theta"))


# Check phi
chagas_summary %>% filter(str_detect(variable, "phi")) %>% pull(median) %>% hist()

chagas_summary %>% filter(str_detect(variable, "E_y"))  %>% arrange(-rhat)
chagas_summary %>% filter(str_detect(variable, "sigma_delta"))

chagas_summary %>% filter(str_detect(variable, "rho"))
chagas_summary %>% filter(str_detect(variable, "alpha"))
chagas_summary %>% filter(str_detect(variable, "delta_pi"))

mcmc_trace(chagas_sample$draws("theta[5,14]"))

plot(chagas_summary$rhat, chagas_summary$ess_bulk)
