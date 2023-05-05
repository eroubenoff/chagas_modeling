# Helper function to create Stan Data -- shared across runtimes

nb_data_funs<- new.env()    
source("nb_data_funs.R", local=nb_data_funs)
attach(nb_data_funs, name="nb_data_funs")

load("chagas_data.Rdata")
# Function to create stan data, which is formatted as a list
# If covariate = FALSE, do not read in climate covariates
# if covariate = TRUE and ncomp = 0, just return the actual variables themselves
# if covariate = TRUE and ncomp >= 1, perform PCA and return the first n prcomps
create_stan_data <- function(
    chagas_arr,  # Array of counts
    pop_arr,     # Array of population
    br_shp,      # Shapefile
    testing=FALSE,  # If TRUE, runs on a subset of the data
    covariate = FALSE, # If TRUE, returns covariates in data object, else returns just the data
    ncomp = 0, # Number of principal components; if 0, returns the covariates  
    pc_obj = FALSE # Return the whole princomp object?
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
  
  if (!testing){
    # There are 3 islands. Need to drop them from the adjacency structure.
    br_nb <- spdep::poly2nb(br_shp)
    # Get them by muni_id to drop from both dfs
    to_drop <- c(1526, 3500)
    br_shp <- br_shp[-to_drop, ]
    chagas_arr <- chagas_arr[, -to_drop]
    pop_arr <- pop_arr[, -to_drop]
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
  
  covariate_arr <- NULL
  nvar <- 0
  
  
  # Add in spatial covariates from the climate files
  if (covariate) {
    covariate_df <- tibble(muni_code = numeric(), year = numeric(), varname = character(), value = numeric())
    
    # Read in temperature data
    path <- "climate_data/CPC/Brazil"
    for (f in list.files(path, full.names = TRUE)) {
      # Calculate the following annual stats: 
      # average tmax, number of days above 25deg, 
      # average tmin
      # average daily precip, number of days over 10mm precip
      
      climate_df <- read_rds(f) %>% select(-area)
      
      # Subset for munis that are actually in our dataset
      climate_df <- climate_df %>% filter(muni_code %in% br_shp$muni_code)
      
      # Extract year from file name
      climate_year <- as.numeric(str_extract(f, "(\\d)+"))
      if (!climate_year %in% yrange) next
      
      # Separate the data 
      muni_code <- climate_df %>% select(muni_code)
      climate_df <- climate_df %>% select(-muni_code)
      
      if (str_detect(f, "max")) {
        # Calculate max stats
        covariate_df <- bind_rows(
          covariate_df,
          bind_cols(
            muni_code = muni_code,
            year = climate_year,
            varname = "avg_tmax",
            value = rowMeans(climate_df, na.rm=TRUE)
          )
        )
        covariate_df <- bind_rows(
          covariate_df,
          bind_cols(
            muni_code = muni_code,
            year = climate_year,
            varname = "n_tmax_gt_30",
            value = rowSums(climate_df > 30, na.rm=TRUE)
          )
        )
        
      }
      else if (str_detect(f, "min")) {
        # Calculate min stats
        covariate_df <- bind_rows(
          covariate_df,
          bind_cols(
            muni_code = muni_code,
            year = climate_year,
            varname = "avg_tmin",
            value = rowMeans(climate_df, na.rm=TRUE)
          )
        )
        
      }
      else if (str_detect(f, "precip")) {
        # Calculate precip stats
        covariate_df <- bind_rows(
          covariate_df,
          bind_cols(
            muni_code = muni_code,
            year = climate_year,
            varname = "avg_precip",
            value = rowMeans(climate_df, na.rm=TRUE)
          )
        )
        
        covariate_df <- bind_rows(
          covariate_df,
          bind_cols(
            muni_code = muni_code,
            year = climate_year,
            varname = "n_precip_gt_10",
            value = rowSums(climate_df > 10, na.rm=TRUE)
          )
        )
      }
    }
    
    # Next pull NDVI
    br_ndvi <- read_csv("br_ndvi.csv")
    br_ndvi <- br_ndvi %>% 
      mutate(date = lubridate::ymd(str_extract(`system:index`, "[\\d+]{4}_[\\d+]{2}_[\\d+]{2}")))
    
    # Calculate summaries
    # Max NDVI, min NDVI, mean NDVI
    br_ndvi <- br_ndvi %>%
      group_by(muni_code, year = lubridate::year(date)) %>% 
      filter(muni_code %in% br_shp$muni_code,
             year %in% yrange) %>% 
      summarize(max_ndvi = max(NDVI, na.rm = TRUE),
                min_ndvi = min(NDVI, na.rm = TRUE),
                mean_ndvi = mean(NDVI, na.rm=TRUE))
    
    
    
    br_ndvi <- br_ndvi %>% 
      pivot_longer(c(max_ndvi, min_ndvi, mean_ndvi), names_to = "varname", values_to = "value")
    

    covariate_df <- bind_rows(covariate_df, br_ndvi)
    
    # For whatever reason there are a few municipalities missing in this dataset
    # There are a handful of incomplete cases (360 of 100k+), and I cannot find why that is
    # the case. Replace those with the overall average.
    
    # Expand to make full set
    covariate_df <- expand_grid(
      year = as.numeric(dimnames(chagas_arr)[[1]]),
      muni_code = as.numeric(dimnames(chagas_arr)[[2]]),
      varname = unique(covariate_df$varname)
    ) %>% 
      left_join(covariate_df)
    
    covariate_df <- covariate_df %>% 
      mutate(value = if_else(varname == "n_precip_gt_10" & value == 0, NA, value),
             value = if_else(varname == "n_precip_gt_10" & value == 0, NA, value))
    
    incompletes <- covariate_df[!complete.cases(covariate_df),]
    covariate_df <- covariate_df[complete.cases(covariate_df),]

    completes <- covariate_df %>% 
      drop_na() %>% 
      group_by(varname)  %>% 
      summarize(value = mean(value))
    
    incompletes <- incompletes %>% 
      select(-value) %>% 
      left_join(completes, by = c("varname"))
    
    covariate_df <- covariate_df %>% 
      bind_rows(incompletes)
    
    
    
    ## Reshape the climate DF into a T,N,variable matrix
    covariate_arr <- covariate_df %>% acast(year ~ muni_code ~ varname, value.var = "value")
    
    if (ncomp >= 1) {
      # Principal components
      covariate_df <- covariate_df %>% 
        pivot_wider(names_from = varname)
        
      covariate_pc <- covariate_df %>% 
        select(-year, -muni_code) %>%
        # drop_na() %>% 
        prcomp(scale=TRUE)
      
      if (pc_obj){
      # If true, return the entire PC object
        covariate_arr <- covariate_pc
      }
      else {
        covariate_pc <- bind_cols(muni_code = covariate_df$muni_code,
                                  year = covariate_df$year,
                                  covariate_pc$x[,1:ncomp])
        
        covariate_arr <- covariate_pc %>%
          pivot_longer(-c(muni_code, year), names_to = "pc", values_to = "value") %>% 
          acast(year ~ muni_code ~ pc, value.var = "value")
        
        ## Make sure all are in the same order as the E and Y arrays 
        covariate_arr <- covariate_arr[dimnames(chagas_arr)[[1]], dimnames(chagas_arr)[[2]], ]
        
        which(!dimnames(covariate_arr)[[2]] %in% dimnames(chagas_arr)[[2]])
        which(!dimnames(chagas_arr)[[2]] %in% dimnames(covariate_arr)[[2]])
        
       nvar <- dim(covariate_arr)[3]
        
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
    nonzero_max = nonzero_max,
    nvar = nvar,
    covariate= covariate_arr
  )
  
  return(stan_data)
}

test_data <- create_stan_data(chagas_arr, pop_arr, br_shp, testing=FALSE, covariate=TRUE, ncomp = 3, pc_obj = FALSE)
test_data <- create_stan_data(chagas_arr, pop_arr, br_shp, testing=FALSE, covariate=TRUE, ncomp = 0)
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

summary(test_data$covariate)
screeplot(test_data$covariate)
biplot(test_data$covariate)
ggplot(aesp)
loadings <- test_data$covariate$rotation %>%
  melt()
colnames(loadings) <- c("varname", "pc", "value")
loadings <- as_tibble(loadings)
ggplot(loadings) + 
  geom_col(aes(x = varname, y = value), fill = if_else(loadings$value > 0, "black", "red")) + 
  facet_wrap(.~pc) + 
  cowplot::theme_cowplot() + 
  coord_flip()

# 
# # Variance explained
# covariate_pc$sdev^2 / sum(covariate_pc$sdev^2)

# Project the first three PCs
# covariate_pc <- predict(covariate_pc, covariate_df, ncomp=3)








