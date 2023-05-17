# Helper function to create Stan Data -- shared across runtimes

nb_data_funs<- new.env()    
source("nb_data_funs.R", local=nb_data_funs)
attach(nb_data_funs, name="nb_data_funs")

# Function to create stan data, which is formatted as a list
# If covariate = FALSE, do not read in climate covariates
# if covariate = TRUE and ncomp = 0, just return the actual variables themselves
# if covariate = TRUE and ncomp >= 1, perform PCA and return the first n prcomps
create_stan_data <- function(
    chagas_arr,  # Array of counts
    pop_arr,     # Array of population
    br_shp,      # Shapefile
    covariate_arr, # if 0, no covariates; else, array of covariates
    testing=FALSE,  # If TRUE, runs on a subset of the data
    nvar= 6 # Number of principal components; if 0, returns the covariates  
    ) {
  
  # If testing, then restrict to just the 3 highest incident states over 
  # 2001- 2005 
  if (testing){
    
    drop_codes <- br_shp %>%
      st_drop_geometry() %>% 
      filter(abbrev_state %in% c("PA", "MA", "AP")) %>% 
      pull(muni_code)
    
    br_shp <- br_shp %>% 
      filter(muni_code %in% drop_codes)
    
    chagas_arr <- chagas_arr[, as.character(drop_codes)]
    pop_arr <- pop_arr[, as.character(drop_codes)]
    
  }

  yrange  <- dimnames(chagas_arr)[[1]]
  if (testing) {
    yrange <- as.character(2001:2005)
    # Limit to first 5 years
    chagas_arr <- chagas_arr[yrange,] 
    pop_arr <- pop_arr[yrange,] 
  }

  # Pull the adjacencies
  br_nb <- spdep::poly2nb(br_shp)
  nbs  <- nb2graph(br_nb)
  N = nbs$N;
  node1 = nbs$node1;
  node2 = nbs$node2;
  N_edges = nbs$N_edges;
  scaling_factor = scale_nb_components(br_nb)[1];
  
  
  ## Tally zero and nonzeros
  n_T = nrow(chagas_arr)
  zero_max = array(rep(0,n_T))
  nonzero_max = array(rep(0,n_T))
  zero_idx = matrix(0, nrow = n_T, ncol = N)
  nonzero_idx = matrix(0, nrow = n_T, ncol = N)
  
  for (t in 1:n_T) {
    zero_max[t] = 0;
    nonzero_max[t] = 0;
    for (n in 1:N){
      if (chagas_arr[t,n] == 0) {
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
    T = length(yrange),
    N_edges = N_edges,
    node1 = node1,
    node2 = node2,
    y = chagas_arr, 
    E = pop_arr, 
    scaling_factor = scaling_factor,
    zero_idx = zero_idx,
    zero_max = zero_max,
    nonzero_idx = nonzero_idx,
    nonzero_max = nonzero_max)
  
  # Add in spatial covariates from the climate 
  if (!missing(covariate_arr)) {
    
    # Extract the corresponding years and PCs dimensions 
    covariate_arr <- covariate_arr[dimnames(chagas_arr)[[1]],dimnames(chagas_arr)[[2]],1:nvar]
    
    stan_data[["nvar"]] <- nvar
    stan_data[["covariate"]] <- covariate_arr
 
  }
 
  
  return(stan_data)
}

# test_data <- create_stan_data(chagas_arr, pop_arr, br_shp, testing=FALSE)
# test_data <- create_stan_data(chagas_arr, pop_arr, br_shp, covariate_arr=covariate_arr, testing=FALSE,  nvar= 3)
# test_data <- create_stan_data(chagas_arr, pop_arr, br_shp, testing=FALSE, covariate=TRUE, ncomp = 0)
# str(test_data$covariate)
# test_data$nvar
# 
# df <- melt(test_data$covariate) 
# colnames(df) <- c("year", "muni_code", "var", "value")
# df <-  df %>% 
#   tibble() %>% 
#   pivot_wider(names_from = var, values_from = value) 
#   # select(-year, -muni_code)  
# df[!complete.cases(df),]

# 

# summary(test_data$covariate)
# screeplot(test_data$covariate)
# biplot(test_data$covariate)
# ggplot(aesp)
# loadings <- test_data$covariate$rotation %>%
#   melt()
# colnames(loadings) <- c("varname", "pc", "value")
# loadings <- as_tibble(loadings)
# ggplot(loadings) + 
#   geom_col(aes(x = varname, y = value), fill = if_else(loadings$value > 0, "black", "red")) + 
#   facet_wrap(.~pc) + 
#   cowplot::theme_cowplot() + 
#   coord_flip()

# 
# # Variance explained
# covariate_pc$sdev^2 / sum(covariate_pc$sdev^2)

# Project the first three PCs
# covariate_pc <- predict(covariate_pc, covariate_df, ncomp=3)








