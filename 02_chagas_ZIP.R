#-------------------------------------------------------------------------------
#### 02_chagas_ZIP.R ####
# Performs a simple ZIP model for Chagas, with no hierarchical parameter
# structure
# Assumes a file "02_chagas_ZIP.stan" with the model code
#-------------------------------------------------------------------------------

library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)

start_time <- Sys.time()

load("chagas_data.Rdata")


stan_data <- list(
  chagas = chagas_count,
  pop = chagas_pop,
  U = dim(chagas_count)[1],
  N_max = dim(chagas_count)[2],
  Y = dim(chagas_count)[3],
  M = dim(chagas_count)[4],
  N = N
)
stan_inits <- function() {
  list(
    p = runif(1, 0, 1),
    l = runif(1, 0, 10),
    lambda = runif(1, 0, 10)
  )
}


chagas_ZIP <- stan("02_chagas_ZIP.stan",
                      model_name = "chagas_ZIP",
                      data = stan_data,
                      init = stan_inits,
                      cores = 4,
                      iter = 1000
                      #, verbose = TRUE
)

save(chagas_ZIP, file = "mcmc_out/chagas_ZIP.Rdata")
end_time <- Sys.time()
message("Model began at ", start_time, " and ended at ", end_time,
        "\nTotal run time is: ", end_time - start_time)


load("mcmc_out/chagas_ZIP.Rdata")
posterior <- mcmc(as.data.frame(chagas_ZIP))
summary(posterior)
mcmc_areas(posterior,
           pars = c("lambda"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("theta"),
           prob = 0.8)

