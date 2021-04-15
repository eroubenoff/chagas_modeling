#-------------------------------------------------------------------------------
#### 04_chagas_offset.R ####
# Spatial offset model (non-covariate)
# Assumes a file "04_chagas_offset.stan" with the model code
# Spatial model only with offset
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html
# https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/bym2_offset_only.stan
# download.file("https://raw.githubusercontent.com/stan-dev/example-models/master/knitr/car-iar-poisson/nb_data_funs.R",
#               destfile= "nb_data_funs.R")
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
source("nb_data_funs.R")



br_nb <- spdep::poly2nb(br_shp)
nbs  <- nb2graph(br_nb)
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;
scaling_factor = scale_nb_components(br_nb)[1];




stan_data <- list(
  N,
  N_edges,
  node1,
  node2,
  y,
  E,
  scaling_factor
)

stan_inits <- function() {
  list(
    p = runif(1, 0, 1),
    l = runif(1, 0, 10),
    lambda = runif(1, 0, 10)
  )
}


chagas_offset <- stan("chagas_offset.stan",
                      model_name = "chagas_offset", 
                      data = stan_data, 
                      control = list(adapt_delta = 0.97), 
                      chains=4, 
                      warmup=7000, 
                      iter=8000, 
                      save_warmup=FALSE)






save(chagas_offset, file = "mcmc_out/chagas_offset.Rdata")
end_time <- Sys.time()


message("Model began at ", start_time, " and ended at ", end_time,
        "\nTotal run time is: ", end_time - start_time)

load("mcmc_out/chagas_offset.Rdata")

posterior <- mcmc(as.data.frame(chagas_offset))

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

