#-------------------------------------------------------------------------------
#### 03_chagas_ZIP_hierarchical.R ####
# Performs a ZIP model with a hierarchical nesting structure---partial pooling
# at the UF (state) level
# Assumes a file "03_chagas_ZIP_hierarchical.stan" with the model code
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
    # p = runif(1, 0, 1),
    # l = runif(1, 0, 10),
    lambda = runif(stan_data$U, 0, 1), 
    theta = runif(stan_data$U, 0, 1)
  )
}

chagas_ZIP_hierarchical <- stan("03_chagas_ZIP_hierarchical.stan",
                      model_name = "chagas_ZIP_hierarchical",
                      data = stan_data,
                      init = stan_inits,
                      cores = 4,
                      iter = 2000
)
save(chagas_ZIP_hierarchical, file = "mcmc_out/chagas_ZIP_hierarchical.Rdata")
end_time <- Sys.time()
message("Model began at ", start_time, " and ended at ", end_time,
        "\nTotal run time is: ", end_time - start_time)
# Chain 1:  Elapsed Time: 21975.6 seconds (Warm-up)
# Chain 1:                25737.5 seconds (Sampling)
# Chain 1:                47713.1 seconds (Total)
  
load("mcmc_out/chagas_ZIP_hierarchical.Rdata")
posterior <- mcmc(as.data.frame(chagas_ZIP_hierarchical))
summary(posterior)
# Plot lambdas against thetas
posterior_df <- data.frame(
  uf_no = 1:27
)
posterior_df$theta <- colMeans(posterior)[1:27]
posterior_df$theta_var <- apply(posterior, 2, var)[1:27]
posterior_df$lambda <- colMeans(posterior)[28:54]
posterior_df$lambda_var <- apply(posterior, 2, var) [28:54]

ggplot(posterior_df) + 
  geom_point(aes(theta, lambda))

mcmc_areas(posterior,
           pars = c("lambda[1]"),
           prob = 0.8)
mcmc_areas(posterior,
           pars = c("theta"),
           prob = 0.8)

library(geobr)

uf_shp <- geobr::read_state()

posterior_shp <- chagas_arr %>% 
  select(uf, uf_no) %>% 
  distinct() %>%
  left_join(posterior_df) %>%
  left_join(br_shp, .)

qtm(posterior_shp, fill = "theta")
qtm(posterior_shp, fill = "lambda")



