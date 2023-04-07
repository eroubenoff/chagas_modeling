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
# we recommend running this is a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
# install_cmdstan(cores = 4)

start_time <- Sys.time()

load("chagas_data.Rdata")
load("./from_ayesha/cleaned_data/pop_all.Rdata")
source("nb_data_funs.R")


br_shp <- br_shp %>% arrange(muni_code)

chagas_offset_count <- chagas_arr %>%
  group_by(year, municipio) %>%
  summarize(count = sum(count))

# Need to create a crosswalk of municipalitity to index
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

# There are 3 islands. Need to drop them from the adjacency structure.
br_nb <- spdep::poly2nb(br_shp)
# Get them by muni_id to drop from both dfs
to_drop <- br_shp %>% slice(1524, 3498, 5564) %>% pull(muni_code)
chagas_offset_count <- chagas_offset_count %>% filter(!municipio %in% to_drop)
br_shp <- br_shp %>%filter(!muni_code %in% to_drop)



# Rectangularize both the pop and chagas counts 
# First make sure both datasets are in the same order
chagas_offset_count <- chagas_offset_count %>% 
  arrange(match(municipio, br_shp$muni_code))

chagas_offset_count <- chagas_offset_count %>% 
  acast(year ~ municipio, value.var = "count") 

chagas_pop <- br_shp %>%
  st_drop_geometry() %>%
  select(muni_code, populacao) %>% 
  slice(rep(1:n(), each = 19)) %>%
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


stan_data <- list(
  N = N,
  T = 19,
  N_edges = N_edges,
  node1 = node1,
  node2 = node2,
  y = chagas_offset_count,
  E = chagas_pop,
  scaling_factor = scaling_factor
)

# stan_inits <- function() {
#   list(
#     p = runif(1, 0, 1),
#     l = runif(1, 0, 10),
#     lambda = runif(1, 0, 10)
#   )
# }
chagas_offset <- cmdstan_model("05_chagas_offset_time.stan",
              cpp_options = list(stan_opencl = TRUE))

chagas_offset$sample(data = stan_data, 
                      # control = list(adapt_delta = 0.99, max_treedepth=20),
                      # pars = c("psi"),
                      chains=4, 
                      parallel_chains=4, 
                      iter_warmup=1000, 
                      iter_sampling=2000, 
                      save_warmup=FALSE)



save(chagas_offset, file = "mcmc_out/chagas_offset.Rdata")
end_time <- Sys.time()


message("Model began at ", start_time, " and ended at ", end_time,
        "\nTotal run time is: ", end_time - start_time)

load("mcmc_out/chagas_offset.Rdata")

posterior <- mcmc(as.data.frame(chagas_offset))


str(posterior)
mcmc_areas(posterior, 
           "psi[1]")

# Extract means to shp
chagas_df$posterior_means = colMeans(posterior)[1:5507]
qtm(chagas_df, fill = "posterior_means")

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

