#-------------------------------------------------------------------------------
#### 04_chagas_offset.R ####
# Spatial offset model (non-covariate), with time
# Assumes a file "05_chagas_offset_time.stan" with the model code
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
# set_cmdstan_path(path="G:/Documents/.cmdstan/cmdstan-2.31.0")

source("nb_data_funs.R")

start_time <- Sys.time()

load("chagas_data.Rdata")
load("./from_ayesha/cleaned_data/pop_all.Rdata")

# Function to create stan data, which is formatted as a list
create_stan_data <- function(chagas_arr, br_shp, testing=FALSE) {
  
  # If testing, then restrict to just the 3 highest incident states over 
  # 2001- 2005 
  if (testing){
    chagas_arr %>% 
      group_by(uf) %>%
      summarize(count = sum(count)) %>%
      arrange(count)
    
    chagas_arr <- chagas_arr %>%
      # filter(uf == "PA") ## PA is highest, RR is lowest
      filter(uf %in% c("PA", "MA", "AP"))
    
    br_shp <- br_shp %>% 
      filter(uf %in% c("PA", "MA", "AP"))
    
  }
  
  br_shp <- br_shp %>% arrange(muni_code)
  
  chagas_offset_count <- chagas_arr %>%
    group_by(year, municipio) %>%
    summarize(count = sum(count))
  
  if (testing) {
    # Limit to first 5 years
    chagas_offset_count <- chagas_offset_count %>%
      filter(year %in% 2001:2005)
  }
  
  # First, confirm that all munis in the chagas DF are in the 
  # shapefile (does not have to be the reverse)
  
  if (!all(
    (chagas_offset_count  %>% pull(municipio))
    %in%
    (br_shp %>% st_drop_geometry() %>% pull(muni_code)) 
  ) ) {
    stop("Not all municipalities in chagas_df found in shapefile")
  }
  
  if (!all(
    (br_shp %>% st_drop_geometry() %>% pull(muni_code)) 
    %in%
    (chagas_offset_count  %>% pull(municipio))
  ) ) {
    stop("Not all municipalities shapefile found in chagas_df")
  }
  
  if (!testing){
    # There are 3 islands. Need to drop them from the adjacency structure.
    br_nb <- spdep::poly2nb(br_shp)
    # Get them by muni_id to drop from both dfs
    to_drop <- br_shp %>% slice(1524, 3498, 5564) %>% pull(muni_code)
    chagas_offset_count <- chagas_offset_count %>% filter(!municipio %in% to_drop)
    br_shp <- br_shp %>%filter(!muni_code %in% to_drop)
  }
  
  
  
  # Rectangularize both the pop and chagas counts 
  # First make sure both datasets are in the same order
  chagas_offset_count <- chagas_offset_count %>% 
    arrange(match(municipio, br_shp$muni_code))
  
  chagas_offset_count <- chagas_offset_count %>% 
    acast(year ~ municipio, value.var = "count") 
  
  chagas_pop <- br_shp %>%
    st_drop_geometry() %>%
    select(muni_code, populacao) %>% 
    slice(rep(1:n(), each = nrow(chagas_offset_count))) %>%
    group_by(muni_code) %>%
    mutate(year = dimnames(chagas_offset_count)[[1]])  %>%
    ungroup()
  
  
  chagas_pop <- chagas_pop %>% 
    mutate(populacao = round(as.numeric(populacao))) %>%
    acast(year ~ muni_code, value.var = "populacao") 
  
  
  # Pull the adjacencies
  br_nb <- spdep::poly2nb(br_shp)
  nbs  <- nb2graph(br_nb)
  N = nbs$N;
  node1 = nbs$node1;
  node2 = nbs$node2;
  N_edges = nbs$N_edges;
  scaling_factor = scale_nb_components(br_nb)[1];
  
  
  ## Tally zero and nonzeros
  n_T = nrow(chagas_offset_count)
  zero_max = array(rep(0,n_T))
  nonzero_max = array(rep(0,n_T))
  zero_idx = matrix(0, nrow = n_T, ncol = N)
  nonzero_idx = matrix(0, nrow = n_T, ncol = N)
  
  for (t in 1:n_T) {
    zero_max[t] = 0;
    nonzero_max[t] = 0;
    for (n in 1:N){
      if (chagas_offset_count[t,n] == 0) {
        zero_max[t] = zero_max[t] +  1;
        zero_idx[t, zero_max[t]] = n;
      }
      else {
        nonzero_max[t] = nonzero_max[t] + 1;
        nonzero_idx[t, nonzero_max[t]] = n;
      }
    }
  }
  
  stan_data <- list(
    N = N,
    T = nrow(chagas_offset_count),
    N_edges = N_edges,
    node1 = node1,
    node2 = node2,
    y = chagas_offset_count, 
    E = chagas_pop, 
    scaling_factor = scaling_factor,
    zero_idx = zero_idx,
    zero_max = zero_max,
    nonzero_idx = nonzero_idx,
    nonzero_max = nonzero_max
  )
  
  return(stan_data)
}

 
# Sys.setenv(STAN_OPENCL=FALSE)
options("cmdstanr_verbose" = TRUE)
# Need to run the following at least once to make sure STAN can use the 
# graphics card
# 
# path_to_opencl_lib <- "G:/CUDA Development/lib/x64"
# cpp_options = list(
#   paste0("LDFLAGS+= -L\"",path_to_opencl_lib,"\" -lOpenCL")
# )
# cmdstanr::cmdstan_make_local(cpp_options = cpp_options)
# cmdstanr::rebuild_cmdstan(cores=4)

# opencl_options = list(
#   stan_opencl = TRUE,
#   opencl_platform_id = 0,
#   opencl_device_id = 0 #in your case its 1 here
# )

testing <- FALSE
stan_data <- create_stan_data(chagas_arr, br_shp, testing = testing)
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
# chagas_sample <- chagas_offset$variational(data = stan_data, 
#                       # chains=ifelse(testing, 2, 4), 
#                       # parallel_chains=4,
#                       # iter_warmup=ifelse(testing, 500, 2000), 
#                       grad_samples = 10,
#                       # elbo_samples = 100,
#                       iter=ifelse(testing,100000, 1000), 
#                       adapt_iter = 1000,
#                       output_dir = "mcmc_out",
#                       save_latent_dynamics = TRUE,
#                       # adapt_delta = 0.8,
#                       refresh=1)


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
