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

setwd("~/Github/chagas_modeling")
library(tidyverse)
library(sf)
library(rstan)
library(StanHeaders)
library(reshape2)
library(bayesplot)
library(coda)
library(tmap)
tmap_mode("view")
# we recommend running this is a fresh R session or restarting your current session
# install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
# install_cmdstan(cores = 4)
set_cmdstan_path(path="G:/Documents/.cmdstan/cmdstan-2.31.0")

source("nb_data_funs.R")

start_time <- Sys.time()

load("chagas_data.Rdata")
load("./from_ayesha/cleaned_data/pop_all.Rdata")

# For testing purposes only:
# The state with the highest count is:
testing <- TRUE 
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


stan_data <- list(
  N = N,
  T = nrow(chagas_offset_count),
  N_edges = N_edges,
  node1 = node1,
  node2 = node2,
  y = chagas_offset_count, 
  E = chagas_pop, 
  scaling_factor = scaling_factor
)

 
Sys.setenv(STAN_OPENCL=TRUE)
options("cmdstanr_verbose" = TRUE)
# Need to run the following at least once to make sure STAN can use the 
# graphics card
# 
# path_to_opencl_lib <- "G:/CUDA Development/lib/x64"
# cpp_options = list(
#   paste0("LDFLAGS+= -L\"",path_to_opencl_lib,"\" -lOpenCL")
# )
# cmdstanr::cmdstan_make_local(cpp_options = cpp_options)
# cmdstanr::rebuild_cmdstan()

chagas_offset <- cmdstan_model("05_chagas_offset_time.stan",
                               # force_recompile = TRUE,
              cpp_options = list(stan_opencl = TRUE))

message("Sampling began at", Sys.time())
chagas_sample <- chagas_offset$sample(data = stan_data, 
                      chains=ifelse(testing, 1, 4), 
                      parallel_chains=ifelse(testing, 1, 4),
                      iter_warmup=ifelse(testing, 500, 2000), 
                      iter_sampling=ifelse(testing,500, 1000), 
                      opencl_ids = c(0, 0),
                      output_dir = "mcmc_out",
                      save_warmup=FALSE)

# chagas_sample$profiles()

chagas_sample$save_output_files("mcmc_out/")
chagas_summary <- chagas_sample$summary()
chagas_summary
# Check convergence
ggplot() + geom_histogram(aes(chagas_summary$rhat)) + geom_vline(aes(xintercept = 1.01))

mean(chagas_summary$ess_bulk)
sd(chagas_summary$ess_bulk)

# Some variables clearly didn't converge:
chagas_summary %>%
  arrange(-rhat) %>% 
  slice_head(prop=0.01)

chagas_summary %>%
  arrange(ess_bulk) %>% 
  slice_head(prop=0.01)

# I think it is mostly phi_lambda when counts are low

chagas_summary %>% filter(str_detect(variable, "psi"))


# Look at some trace plots
mcmc_areas(chagas_sample$draws("beta0"))
mcmc_areas(chagas_sample$draws("alpha0"))
mcmc_areas(chagas_sample$draws("rho_pi"))
mcmc_areas(chagas_sample$draws("rho_lambda"))
mcmc_areas(chagas_sample$draws("sigma_delta"))


chagas_sample$save_object(file="mcmc_out/chagas_offset.RDS")
chagas_sample <- readRDS("mcmc_out/chagas_offset.RDS")
#save(chagas_sample, file = "mcmc_out/chagas_offset.Rdata")
end_time <- Sys.time()

# Pick a few municipalities and plot their parameters over time
chagas_summary %>% 
  filter(str_detect(variable, "^pi\\[[0-9]*,12\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) + geom_errorbar(aes(x = t,ymin = q5, ymax = q95))

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) + geom_errorbar(aes(x = t,ymin = q5, ymax = q95))


# Test mapping
IR <- chagas_summary %>%
  filter(str_detect(variable, "^E_y\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")
IR <- IR %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

count <- as.data.frame(t(chagas_offset_count[1:5, ])) 
colnames(count) = c("MLE_y1","MLE_y2", "MLE_y3", "MLE_y4", "MLE_y5")

IR <- bind_cols(IR, count) 
IR <- bind_cols(br_shp, IR)

IR <- IR %>% 
  mutate(populacao = as.numeric(populacao)) %>%
  mutate(across(c(year1, year2, year3, year4, year4, MLE_y1, MLE_y2, MLE_y3, MLE_y4, MLE_y5), ~./populacao))


tmap_arrange(
  tm_shape(IR) + tm_polygons(col = "year5", style = "cont"),
  tm_shape(IR) + tm_polygons(col = "MLE_y5", style = "cont")
)






pi <- chagas_summary %>%
  filter(str_detect(variable, "^pi\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")

pi <- pi %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(x))
}
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

message("Model began at ", start_time, " and ended at ", end_time,
        "\nTotal run time is: ", end_time - start_time)


