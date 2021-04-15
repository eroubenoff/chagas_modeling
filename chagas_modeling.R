# Chagas Modeling 
library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)

start_time <- Sys.time()

#-------------------------------------------------------------------------------
#### Load Data ####
#-------------------------------------------------------------------------------

setwd("~/90days/eroubenoff/chagas")
load("./chagas_data/chagas_municipio.Rdata")
br_shp <- st_read("./from_ayesha/municipios_2010.shp", "municipios_2010")
load("./from_ayesha/cleaned_data/pop_all.Rdata")

#-------------------------------------------------------------------------------
#### Extract code from shapefile ####
#-------------------------------------------------------------------------------
br_shp <- br_shp %>%
  mutate(muni_code = as.numeric(substr(codigo_ibg, 1, 6)))





#-------------------------------------------------------------------------------
#### Rectangularize the Chagas data ####
#-------------------------------------------------------------------------------
# Pivot and extract muni code
chagas_municipio <- chagas_municipio %>%
  pivot_longer(cols = 2:14, names_to = "month", values_to = "count") %>%
  mutate(muni_code = str_sub(`Município de residência`, 1, 6),
         muni_code = as.numeric(muni_code))

# Create df
chagas_arr <- expand_grid(
  municipio = unique(br_shp$muni_code),
  year = unique(chagas_municipio$year),
  month = unique(chagas_municipio$month)
)

# Join df with data
chagas_arr <- left_join(chagas_arr,
                        chagas_municipio %>% select(
                          municipio = muni_code,
                          year, month, count),
)

# Add in UF 
chagas_arr <- left_join(chagas_arr,
                        br_shp %>%
                          st_drop_geometry() %>%
                          select(municipio = muni_code, uf))

# Set NAs to 0
chagas_arr[is.na(chagas_arr)] <- 0

# Lookup indices for municipios and UFs
chagas_arr <- left_join(
  chagas_arr,
  chagas_arr %>% 
    select(uf, municipio) %>%
    distinct() %>%
    arrange(uf, municipio) %>%
    group_by(uf) %>%
    mutate(uf_no = cur_group_id(),
           muni_no = 1:n())
) 
 
  
  
# Cast to array
# Need to order the month column like: factor(ColumnName, levels = unique(ColumnName)
chagas_count <- acast(chagas_arr, formula = 
                        uf_no~muni_no~year~factor(month, levels = unique(month)), 
                      value.var = "count",
                      fill = 0)
# dim(chagas_count)
# dimnames(chagas_count)

# Vector of the number of municipalities within each UF
N <- chagas_arr %>%
  select(uf_no, muni_no) %>%
  distinct() %>%
  group_by(uf_no) %>%
  slice_max(muni_no, n = 1) %>%
  pull(muni_no)
 




#-------------------------------------------------------------------------------
#### Rectangularize the Pop data ####
#-------------------------------------------------------------------------------
# Create pop df
chagas_pop <- chagas_arr %>%
  select(municipio, year, uf_no, muni_no) %>%
  distinct()

# Join df with data
chagas_pop <- chagas_pop %>%
  left_join(pop.all %>% select(municipio = code, year, ps))

# Cast matrix
chagas_pop <- acast(chagas_pop, formula = uf_no~muni_no~year, value.var = "ps", fill = 0)
# dim(chagas_pop)
# dimnames(chagas_pop)

# The population data only go through 2015 so for now we need to limit
# the incidence data to 2015 (NTS: get the rest of the pop data)
chagas_count <- chagas_count[, , 0:15, 1:12]
chagas_pop <- chagas_pop[, ,1:15] 

# See which data are missing by lining up the pop and the count data
if (FALSE) {
  for (u in 1:dim(chagas_pop)[1]) {
    message(u, "\n")
    for (m in 1:dim(chagas_pop)[2]) {
      for (y in 1:dim(chagas_pop)[3]) {
        if (is.na(chagas_pop[u, m, y]) & all(is.na(chagas_count[u, m, y, ]))){
          next
        }
        if (!is.na(chagas_pop[u, m, y]) & all(!is.na(chagas_count[u, m, y, ]))){
          next
        }
        else {
          message("Mismatch at: UF ", u, ", M ", m, ", Y ", y, "missmatch \n")
        }
      }
    }
  }
}

# There are a few missing populations
chagas_pop[12, 35, c(1,2,3)] <- rep(chagas_pop[12, 35, 4], 3)
chagas_pop[13, 51, c(1,2,3)] <- rep(chagas_pop[13, 51, 4], 3)
chagas_pop[13, 52, c(1,2,3)] <- rep(chagas_pop[13, 52, 4], 3)
chagas_pop[17, 14, c(1,2,3)] <- rep(chagas_pop[17, 14, 4], 3)



#-------------------------------------------------------------------------------
#### Modeling ####
#-------------------------------------------------------------------------------
# Example: https://mc-stan.org/docs/2_20/stan-users-guide/zero-inflated-section.html
stan_code <- "
data {
  int<lower=0> U;                       // Number of UFs
  int<lower=0> N_max;                   // Maximum number of municipios
  int<lower=0> N[U];                    // Number of municipios within UF u 
  int<lower=0> Y;                       // Year
  int<lower=0> M;                       // Month
  int<lower=0> chagas[U, N_max, Y, M];  // Array of counts
  int<lower=0> pop[U, N_max, Y];        // Array of Population
}
parameters {
  real<lower=0, upper=1> theta;     // Probability of 0 count
  real<lower=0> lambda;             // Poisson mean P(chagas count | theta == 0)
}
model {
  theta ~ uniform(0, 1);            // Bounded by 0 and 1
  lambda ~ uniform(0, 100);         // Greater than 0
  
  for (u in 1:U) {
    for (i in 1:N[u]) {
      for (y in 1:Y) {
        for (m in 1:M) {
          if (chagas[u, i, y, m] == 0) 
                  target += log_sum_exp(bernoulli_lpmf(1 | theta),
                                  bernoulli_lpmf(0 | theta)
                                + poisson_lpmf(chagas[u, i, y, m] | lambda * pop[u, i, y]));
          else
                  target += bernoulli_lpmf(0 | theta)
                             + poisson_lpmf(chagas[u, i, y, m] | lambda * pop[u, i, y]);             
        }
      }
    }
  }
}
"

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


if (FALSE) {
  chagas_model1 <- stan(model_name = "chagas_model1",
       model_code = stan_code,
       data = stan_data,
       init = stan_inits,
       cores = 4,
       iter = 1000
       #, verbose = TRUE
       )
  save(chagas_model1, file = "mcmc_out/chagas_model1.Rdata")
  end_time <- Sys.time()
  print(end_time - start_time)
}
if(FALSE) {
  load("mcmc_out/chagas_model1.Rdata")
  posterior <- mcmc(as.data.frame(chagas_model1))
  summary(posterior)
  mcmc_areas(posterior,
             pars = c("lambda"),
             prob = 0.8)
  mcmc_areas(posterior,
             pars = c("theta"),
             prob = 0.8)

}





# Model 2: 
# Nesting of municipios within UFs
stan_code <- "
data {
  int<lower=0> U;                       // Number of UFs
  int<lower=0> N_max;                   // Maximum number of municipios
  int<lower=0> N[U];                    // Number of municipios within UF u 
  int<lower=0> Y;                       // Year
  int<lower=0> M;                       // Month
  int<lower=0> chagas[U, N_max, Y, M];  // Array of counts
  int<lower=0> pop[U, N_max, Y];        // Array of Population
}
parameters {
  vector<lower=0, upper=1>[U] theta;     // Probability of 0 count
  vector<lower=0>[U] lambda;             // Poisson mean P(chagas count | theta == 0)
}
model {
  theta ~ uniform(0, 1);            // Bounded by 0 and 1
  lambda ~ uniform(0, 100);         // Greater than 0
  
  for (u in 1:U) {
    // Lambda and theta within each UF
    // theta[u] ~ uniform(0, 1);
    // lambda[u] ~ uniform(0, 100);
    
    for (i in 1:N[u]) {
      for (y in 1:Y) {
        for (m in 1:M) {
          if (chagas[u, i, y, m] == 0) 
                  target += log_sum_exp(bernoulli_lpmf(1 | theta[u]),
                                  bernoulli_lpmf(0 | theta[u])
                                + poisson_lpmf(chagas[u, i, y, m] | lambda[u] * pop[u, i, y]));
          else
                  target += bernoulli_lpmf(0 | theta[u])
                             + poisson_lpmf(chagas[u, i, y, m] | lambda[u] * pop[u, i, y]);             
        }
      }
    }
  }
}
"


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

if (TRUE) {
  chagas_model2 <- stan(model_name = "chagas_model2",
                        model_code = stan_code,
                        data = stan_data,
                        init = stan_inits,
                        cores = 4,
                        iter = 2000
  )
  save(chagas_model2, file = "mcmc_out/chagas_model2.Rdata")
  end_time <- Sys.time()
  print(end_time - start_time)
  # Chain 1:  Elapsed Time: 21975.6 seconds (Warm-up)
  # Chain 1:                25737.5 seconds (Sampling)
  # Chain 1:                47713.1 seconds (Total)
}
if(FALSE) {
  load("mcmc_out/chagas_model2.Rdata")
  posterior <- mcmc(as.data.frame(chagas_model2))
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
  
  
}








# Spatial model only with offset
# https://mc-stan.org/users/documentation/case-studies/icar_stan.html
# https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/bym2_offset_only.stan

download.file("https://raw.githubusercontent.com/stan-dev/example-models/master/knitr/car-iar-poisson/nb_data_funs.R",
              destfile= "nb_data_funs.R")
source("nb_data_funs.R")



br_nb <- spdep::poly2nb(br_shp)
nbs  <- nb2graph(br_nb)
N = nbs$N;
node1 = nbs$node1;
node2 = nbs$node2;
N_edges = nbs$N_edges;
scaling_factor = scale_nb_components(br_nb)[1];

stan_code <-  "
functions {
  real icar_normal_lpdf(vector phi, int N, int[] node1, int[] node2) {
    return -0.5 * dot_self(phi[node1] - phi[node2]);
 }
}
data {
  int<lower=0> N;
  int<lower=0> N_edges;
  int<lower=1, upper=N> node1[N_edges];  // node1[i] adjacent to node2[i]
  int<lower=1, upper=N> node2[N_edges];  // and node1[i] < node2[i]

  int<lower=0> y[N];              // count outcomes
  vector<lower=0>[N] x;           // coefficient
  vector<lower=0>[N] E;           // exposure

}
transformed data {
  vector[N] log_E = log(E);
}
parameters {
  real beta0;             // intercept
  real beta1;             // slope
  real<lower=0> sigma;    // overall standard deviation
  vector[N] phi;         // spatial effects
}
model {
  y ~ poisson_log(log_E + beta0 + beta1 * x + phi * sigma);
  beta0 ~ normal(0.0, 1.0);
  beta1 ~ normal(0.0, 1.0);
  sigma ~ normal(0.0, 1.0);
  phi ~ icar_normal_lpdf(N, node1, node2);
  // soft sum-to-zero constraint on phi
  // more efficient than mean(phi) ~ normal(0, 0.001)
  sum(phi) ~ normal(0, 0.001 * N);
}
generated quantities {
  vector[N] eta = log_E + beta0 + beta1 * x + phi * sigma;
  vector[N] mu = exp(eta);
}
"

chagas_model3 <- stan(model_name = "chagas_model3", 
                      model_code = chagas_model3,
                      data = list(
                        N,
                        N_edges,
                        node1,
                        node2,
                        y,
                        E,
                        scaling_factor
                      ), 
                      control = list(adapt_delta = 0.97), 
                      chains=4, 
                      warmup=7000, 
                      iter=8000, 
                      save_warmup=FALSE)
















