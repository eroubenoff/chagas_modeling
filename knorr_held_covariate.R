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



load("chagas_data.Rdata")
load("./from_ayesha/cleaned_data/pop_all.Rdata")
options("cmdstanr_verbose" = TRUE)

source("chagas_helper.R")
testing <- FALSE 
stan_data <- create_stan_data(chagas_arr, br_shp, testing = testing, covariate=TRUE, ncomp = 3)
chagas_offset <- cmdstan_model("knorr_held_covariate.stan",
                               force_recompile = FALSE)

message("Sampling began at", Sys.time())
chagas_sample <- chagas_offset$sample(data = stan_data,
                                      chains=ifelse(testing, 2, 4),
                                      parallel_chains=4,
                                      iter_warmup=ifelse(testing, 500, 1000),
                                      iter_sampling=ifelse(testing,100, 500),
                                      output_dir = "mcmc_out",
                                      save_latent_dynamics = TRUE,
                                      adapt_delta = 0.8,
                                      refresh=25,
                                      save_warmup=FALSE)

# chagas_summmary <- chagas_sample$summary()
# chagas_summmary
# 
# chagas_summmary %>% filter(str_detect(variable, "beta"))
