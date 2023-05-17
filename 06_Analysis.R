# Chagas Analyis
setwd("~/chagas_modeling")
library(tidyverse)
library(sf)
library(reshape2)
library(tmap)
tmap_mode("plot")
library(posterior)
library(cmdstanr)
library(geobr)
library(cowplot)
library(foreach)
library(doParallel)
# ndvi <- read_csv("~/Downloads/test_ndvi.csv")


load("chagas_data.Rdata")


# Get map of UFs
ufs <- geobr::read_state()

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(-x))
}

## Skip for now
if (FALSE) {
    
  
  # Map of crude case counts
  counts_map <- bind_cols(br_shp, count = chagas_arr %>% colSums() %>% na_if(0)) %>% 
    tm_shape() + 
    tm_fill(col = "count", style= "cont")  +
    tm_layout(aes.palette = list(seq = "-RdYlGn"))  +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'Cumulative Cases between 2001 and 2009', 
              main.title.size = 1)
  
  tmap_save(counts_map, "counts_map.png", height = 5, width = 5)
  
  # Quick analysis: goodness-of-fit 
  chagas_count_v <- as.vector(chagas_offset_count)
  pois_ll <- function(p) {
    sum(-dpois(chagas_count_v, p, log=TRUE))
  }
  pois_ll <- optim(c(0.01), pois_ll, method = "L-BFGS-B", lower = 2e-16, upper = 10)
  
  zip_ll <- function(x) {
    pi <- x[1]
    lambda <- x[2]
    ll = 0
    
    for (i in chagas_count_v){
      if (i == 0){
        ll = ll + -log(pi + (1-pi)*dpois(0,lambda))
      }
      else {
        ll = ll + -log(1-pi) + -(dpois(i, lambda, log=TRUE))
      }
    }
    return(ll)
  }
  zip_ll <- optim(c(0.9,0.01), zip_ll, method = "L-BFGS-B", lower = c(0, 0.0000001), upper = c(0.99999, 10))
  
  lr <- -2*(zip_ll$value - pois_ll$value)
  dchisq(lr, length(chagas_count_v) - 1)
  v <- 1:200000
  ggplot() + 
    geom_line(aes(v, dchisq(v,length(chagas_count_v) - 1 ))) + 
    geom_vline(aes(xintercept =  lr))
  
  
  # Analysis of time trend
  chagas_offset_count %>% group_by(year) %>% summarize(count= sum(count)) %>%
    ggplot() + geom_point(aes(x=year, y=count))
}

nCores <- 4
cluster <- makeCluster(nCores, type = "FORK")
registerDoParallel(cluster)
process_mcmc <- function(fname, out_prefix, other_vars = c(), load = TRUE, draws = FALSE) {
    
  f1 <- str_replace(fname, "X", "1")
  f2 <- str_replace(fname, "X", "2")
  f3 <- str_replace(fname, "X", "3")
  f4 <- str_replace(fname, "X", "4")
  
  # Predicted counts
  pi_summary_f <- paste0(out_prefix, "pi_summary.RDS")
  pi_draws_f <- paste0(out_prefix, "pi_draws.RDS")
  lambda_summary_f <- paste0(out_prefix, "lambda_summary.RDS")
  lambda_draws_f <- paste0(out_prefix, "lambda_draws.RDS")
  delta_lambda_summary_f <- paste0(out_prefix, "delta_lambda_summary.RDS")
  delta_lambda_draws_f <- paste0(out_prefix, "delta_lambda_draws.RDS")
  delta_pi_summary_f <- paste0(out_prefix, "delta_pi_summary.RDS")
  delta_pi_draws_f <- paste0(out_prefix, "delta_pi_draws.RDS")
  phi_alpha_summary_f <- paste0(out_prefix, "phi_alpha_summary.RDS")
  phi_alpha_draws_f <- paste0(out_prefix, "phi_alpha_draws.RDS")
  
  if (!load) {
    
    ## ------
    ## Pi: Load csvs
    ## ------
    message(Sys.time(), "Loading pi1")
    pi1 <- read_cmdstan_csv(f1,variables = "pi", format = "draws_list")
    message(Sys.time(), "Loading pi2")
    pi2 <- read_cmdstan_csv(f2,variables = "pi", format = "draws_list")
    message(Sys.time(), "Loading pi3")
    pi3 <- read_cmdstan_csv(f3, variables = "pi", format = "draws_list")
    message(Sys.time(), "Loading pi4")
    pi4 <- read_cmdstan_csv(f4, variables = "pi", format = "draws_list")
    
    ## ------
    ## Pi: Bind draws 
    ## ------
    message(Sys.time(), "Binding pi")
    pi_draws <- bind_draws(x = pi1$post_warmup_draws, pi2$post_warmup_draws, 
                      pi3$post_warmup_draws, pi4$post_warmup_draws, along = "chain")
    rm(pi1, pi2, pi3, pi4)
    gc()
    
    
    ## ------
    ## Pi: Extract summary 
    ## ------
    message(Sys.time(), "Calculating pi summary")
    pi_names <- variables(pi_draws)
    pi_summary <- foreach(x = pi_names) %dopar% {
      d <- cbind(pi_draws[[1]][[x]], pi_draws[[2]][[x]], pi_draws[[3]][[x]], pi_draws[[4]][[x]])
      c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
    }
    pi_summary <- bind_rows(pi_summary)
    
    ## ------
    ## Pi: Save objs 
    ## ------
    saveRDS(pi_summary, file=pi_summary_f)
    saveRDS(pi_draws, file=pi_draws_f, compress=TRUE)
    rm(pi_draws, pi_summary)
    gc()
    
    
    ## ------
    ## Lambda : Load csvs
    ## ------
    message(Sys.time(), "Loading lambda1")
    lambda1 <- read_cmdstan_csv(f1,
      variables = "lambda", format = "draws_list")
    message(Sys.time(), "Loading lambda2")
    lambda2 <- read_cmdstan_csv(f2,
      variables = "lambda", format = "draws_list")
    message(Sys.time(), "Loading lambda3")
    lambda3 <- read_cmdstan_csv(f3,
      variables = "lambda", format = "draws_list")
    message(Sys.time(), "Loading lambda4")
    lambda4 <- read_cmdstan_csv(f4,
      variables = "lambda", format = "draws_list")
    
    ## ------
    ## Lambda: Bind draws 
    ## ------
    message(Sys.time(), "Binding lambda")
    lambda_draws <- bind_draws(x = lambda1$post_warmup_draws, lambda2$post_warmup_draws, 
                     lambda3$post_warmup_draws, lambda4$post_warmup_draws, along = "chain")
    rm(lambda1, lambda2, lambda3, lambda4)
    gc()
    
    ## ------
    ## Lambda: Extract summary 
    ## ------
    lambda_names <- variables(lambda_draws)
    message(Sys.time(), "Lambda summary")
    lambda_summary <- foreach(x = lambda_names) %dopar% {
      d <- cbind(lambda_draws[[1]][[x]], lambda_draws[[2]][[x]], lambda_draws[[3]][[x]], lambda_draws[[4]][[x]])
      c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
    }
    lambda_summary <- bind_rows(lambda_summary)
    
    ## ------
    ## Lambda: Save objs 
    ## ------
    saveRDS(lambda_draws, lambda_draws_f)
    saveRDS(lambda_summary, lambda_summary_f)
    rm(lambda_summary)
    rm(lambda_draws)
    gc()
    
    
    ## ------
    ## delta : Load csvs
    ## ------
    message(Sys.time(), "Loading delta1")
    delta1 <- read_cmdstan_csv(f1,
                               variables = c("delta_lambda"), format = "draws_list")
    message(Sys.time(), "Loading delta2")
    delta2 <- read_cmdstan_csv(f2,
                               variables = c("delta_lambda"), format = "draws_list")
    message(Sys.time(), "Loading delta3")
    delta3 <- read_cmdstan_csv(f3,
                               variables = c("delta_lambda"), format = "draws_list")
    message(Sys.time(), "Loading delta4")
    delta4 <- read_cmdstan_csv(f4,
                               variables = c("delta_lambda"), format = "draws_list")
    
    ## ------
    ## delta: Bind draws 
    ## ------
    message(Sys.time(), "Binding delta")
    delta_draws <- bind_draws(x = delta1$post_warmup_draws, delta2$post_warmup_draws, 
                              delta3$post_warmup_draws, delta4$post_warmup_draws, along = "chain")
    rm(delta1, delta2, delta3, delta4)
    gc()
    
    ## ------
    ## delta: Extract summary 
    ## ------
    delta_names <- variables(delta_draws)
    message(Sys.time(), "delta summary")
    delta_summary <- foreach(x = delta_names) %dopar% {
      d <- cbind(delta_draws[[1]][[x]], delta_draws[[2]][[x]], delta_draws[[3]][[x]], delta_draws[[4]][[x]])
      c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
    }
    delta_summary <- bind_rows(delta_summary)
    
    ## ------
    ## delta: Save objs 
    ## ------
    saveRDS(delta_draws, delta_lambda_draws_f)
    saveRDS(delta_summary, delta_lambda_summary_f)
    rm(delta_summary)
    rm(delta_draws)
    gc()
    
    
    
    ## ------
    ## delta : Load csvs
    ## ------
    message(Sys.time(), "Loading delta1")
    delta1 <- read_cmdstan_csv(f1,
      variables = c("delta_pi"), format = "draws_list")
    message(Sys.time(), "Loading delta2")
    delta2 <- read_cmdstan_csv(f2,
      variables = c("delta_pi"), format = "draws_list")
    message(Sys.time(), "Loading delta3")
    delta3 <- read_cmdstan_csv(f3,
      variables = c("delta_pi"), format = "draws_list")
    message(Sys.time(), "Loading delta4")
    delta4 <- read_cmdstan_csv(f4,
      variables = c("delta_pi"), format = "draws_list")
    
    ## ------
    ## delta: Bind draws 
    ## ------
    message(Sys.time(), "Binding delta")
    delta_draws <- bind_draws(x = delta1$post_warmup_draws, delta2$post_warmup_draws, 
                     delta3$post_warmup_draws, delta4$post_warmup_draws, along = "chain")
    rm(delta1, delta2, delta3, delta4)
    gc()
    
    ## ------
    ## delta: Extract summary 
    ## ------
    delta_names <- variables(delta_draws)
    message(Sys.time(), "delta summary")
    delta_summary <- foreach(x = delta_names) %dopar% {
      d <- cbind(delta_draws[[1]][[x]], delta_draws[[2]][[x]], delta_draws[[3]][[x]], delta_draws[[4]][[x]])
      c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
    }
    delta_summary <- bind_rows(delta_summary)
    
    ## ------
    ## delta: Save objs 
    ## ------
    saveRDS(delta_draws, delta_pi_draws_f)
    saveRDS(delta_summary, delta_pi_summary_f)
    rm(delta_summary)
    rm(delta_draws)
    gc()
    
    
    ## ------
    ## Others: Load csvs
    ## phi, alpha, theta, mu
    ## ------
    
    other_vars <- c(
      "mu_pi",
      "mu_lambda",
      "alpha_pi",
      "alpha_lambda",
      "sigma_alpha_pi",
      "sigma_alpha_lambda",
      "phi_pi",
      # "phi_lambda",
      "theta_pi",
      "theta_lambda",
      "rho_pi",
      # "rho_lambda",
      "sigma_convolved_pi",
      "sigma_convolved_lambda",
      "sigma_delta_lambda",
      "sigma_delta_pi",
      other_vars
    )
    # Repeat for Phi, alpha, and theta, which are much samller 
    
    message(Sys.time(), "Loading phi_alpha1")
    phi_alpha1 <- read_cmdstan_csv(f1,
      variables = other_vars, format = "draws_list")
    message(Sys.time(), "Loading phi_alpha2")
    phi_alpha2 <- read_cmdstan_csv(f2,
      variables =  other_vars, format = "draws_list")
    message(Sys.time(), "Loading phi_alpha3")
    phi_alpha3 <- read_cmdstan_csv(f3,
      variables =  other_vars, format = "draws_list")
    message(Sys.time(), "Loading phi_alpha4")
    phi_alpha4 <- read_cmdstan_csv(f4,
      variables =  other_vars, format = "draws_list")
    
    
    ## ------
    ## Others: Bind draws 
    ## ------
    message(Sys.time(), "Binding phi alpha")
    phi_alpha_draws <- bind_draws(x = phi_alpha1$post_warmup_draws, phi_alpha2$post_warmup_draws, 
                      phi_alpha3$post_warmup_draws, phi_alpha4$post_warmup_draws, along = "chain")
    
    rm(phi_alpha1, phi_alpha2, phi_alpha3, phi_alpha4)
    gc()
    
    ## ------
    ## Others: Calculate Summary 
    ## ------
    phi_alpha_names <- variables(phi_alpha_draws)
    message(Sys.time(), "phi_alpha summary")
    phi_alpha_summary <- foreach(x = phi_alpha_names) %dopar% { 
      d <- cbind(phi_alpha_draws[[1]][[x]], phi_alpha_draws[[2]][[x]], phi_alpha_draws[[3]][[x]], phi_alpha_draws[[4]][[x]])
      c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
    }
    phi_alpha_summary <- bind_rows(phi_alpha_summary)

    
    ## ------
    ## Others: Save objs 
    ## ------
    saveRDS(phi_alpha_draws, phi_alpha_draws_f)
    saveRDS(phi_alpha_summary, phi_alpha_summary_f)
    
    rm(phi_alpha_draws, phi_alpha_summary)
    
  
  }
  
  ret <- list() 
  ret[["pi_summary"]] <- readRDS(pi_summary_f)
  ret[["lambda_summary"]] <- readRDS(lambda_summary_f)
  ret[["delta_lambda_summary"]] <- readRDS(delta_lambda_summary_f)
  ret[["delta_pi_summary"]] <- readRDS(delta_pi_summary_f)
  ret[["phi_alpha_summary"]] <- readRDS(phi_alpha_summary_f)
  
  if (draws) {
    ret[["pi_draws"]] <- readRDS(pi_draws_f)
    ret[["lambda_draws"]] <- readRDS(lambda_draws_f)
    ret[["delta_lambda_draws"]] <- readRDS(delta_lambda_draws_f)
    ret[["delta_pi_draws"]] <- readRDS(delta_pi_draws_f)
    ret[["phi_alpha_draws"]] <- readRDS(phi_alpha_draws_f)
    
  }
  
  return(ret)
}



load <- FALSE 

fname <- "mcmc_out/knorr_held_convolved_both_parts-202305111148-X-64af77.csv"

main_mcmc <- process_mcmc(fname = fname, out_prefix = "", load = load)


fname <- "mcmc_out/knorr_held_covariate-202305091557-X-61b29d.csv"
covariate_mcmc <- process_mcmc(fname = fname, out_prefix = "covariate_", 
                               other_vars = c("beta_pi", "beta_lambda", 
                                              "sigma_beta_pi", "sigma_beta_lambda"), 
                               load = load)

stopCluster(cluster)
# attach(covariate_mcmc)


pi_summary <- main_mcmc[["pi_summary"]]
lambda_summary <- main_mcmc[["lambda_summary"]]
phi_alpha_summary <- main_mcmc[["phi_alpha_summary"]]


pi_summary <- pi_summary %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,",".")) %>%
  separate(variable, sep = "\\.", into = c("variable","year", "muni", NA))

lambda_summary <- lambda_summary %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,",".")) %>%
  separate(variable, sep = "\\.", into = c("variable","year", "muni", NA))

phi_alpha_summary <- phi_alpha_summary %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) 

muni_crosswalk <- br_shp %>%
  st_drop_geometry() %>%
  select(name_muni, uf = abbrev_state, code_state, muni_code, contains("pop")) %>%
  # mutate(populacao = as.numeric(populacao)) %>%
  rownames_to_column(var = "muni")

muni_crosswalk <- muni_crosswalk %>% as_tibble() %>% 
  pivot_longer(-c(muni, name_muni, uf, code_state, muni_code), values_to = "pop") %>%
  separate(name, into = c(NA, "year"), convert=TRUE)


# Generate E_y from pi and Lambda
E_y  <- pi_summary  %>%
  mutate(variable = str_replace(variable, "pi", "E_y")) %>%
  select(variable, year, muni, median) %>% 
  mutate(year = as.numeric(year) + 2000)
E_y$median  <- (1-pi_summary$median)*lambda_summary$median
E_y
E_y <- E_y %>% left_join(muni_crosswalk, by = c("muni", "year"))
E_y


## Calculate RMSE
main_rmse <- left_join(
  E_y,
  br_shp %>%
    st_drop_geometry() %>% 
    select(muni_code, starts_with("chagas")) %>%
    pivot_longer(-muni_code, values_to = "chagas") %>% 
    separate(name, into = c(NA, "year"), convert=TRUE) 
)

main_rmse %>% 
  # group_by(year) %>%
  summarize(rmse = sqrt(mean((median - chagas)^2)))

##Calculate incidence rate as sum(E_y)/sum(pop)

total_IR <- E_y %>% 
  summarize(IR = 100000*sum(median)/sum(pop))
total_IR 

municipality_IR <- E_y %>%
  group_by(muni_code, uf, name_muni) %>%
  summarize(IR = sum(median)/sum(pop))

municipality_IR

UF_IR <- E_y %>%
  group_by(uf) %>%
  summarize(IR = sum(median)/sum(pop))

UF_IR

# States with the highest IR
UF_IR <- UF_IR %>%
  mutate(log_IR = log(IR),
         IR_per_100k = IR * 100000) %>%
  arrange(-IR)

UF_IR

ufs <- ufs %>%
  left_join(UF_IR, by = c("abbrev_state" = "uf"))

UF_IR_map <- tm_shape(ufs) + 
  tm_polygons(col = "log_IR", style = "cont") +
  tm_text("abbrev_state", size = 1, col = "black") +
  tm_layout(aes.palette = list(seq = "-RdYlGn"))  
  

tmap_save(UF_IR_map, "UF_IR_map.png")

## Look at Para and Amapa 

PA_AP <- E_y %>%
  filter(uf %in% c("AP", "PA")) %>% 
  mutate(year = as.numeric(year))  

# IR
PA_AP %>% 
  group_by(name_muni) %>% 
  summarize(IR = 100000*sum(median)/sum(pop),
            pop = mean(pop)) %>%
  arrange(-IR) %>% View

PA_AP %>% 
  group_by(name_muni) %>% 
  summarize(IR = sum(median),
            pop = mean(pop))%>%
  arrange(-IR)


## Time trend: alpha parameters
alpha <- phi_alpha_summary %>%
  filter(str_detect(variable, "^alpha_[pi|lambda]")) %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,",".")) %>%
  separate(variable, sep = "\\.", into = c("variable","year", NA)) %>% 
  mutate(across(c(year, median, rhat, ess), ~as.numeric(.)))

# Need to add in sigma and mu and convert to the appropriate scale
mu_lambda <- phi_alpha_summary %>%
  filter(str_detect(variable, "mu_lambda"))
mu_lambda <- as.numeric(mu_lambda %>% pull(median))
mu_lambda

mu_pi <- phi_alpha_summary %>%
  filter(str_detect(variable, "mu_pi"))
mu_pi <- as.numeric(mu_pi %>% pull(median))
mu_pi

sigma_alpha_pi <- phi_alpha_summary %>%
  filter(variable == "sigma_alpha_pi")
sigma_alpha_pi <- as.numeric(sigma_alpha_pi %>% pull(median))
sigma_alpha_pi

sigma_alpha_lambda <- phi_alpha_summary %>%
  filter(variable == "sigma_alpha_lambda")
sigma_alpha_lambda <- as.numeric(sigma_alpha_lambda %>% pull(median))
sigma_alpha_lambda

alpha <- alpha %>% 
  mutate(median_transformed = case_when(
    variable == "alpha_pi" ~  sigmoid(mu_pi +median),
    variable == "alpha_lambda"  ~ exp(mu_lambda +  median)
  ))

p1 <- alpha %>%
  ggplot() + 
  geom_line(aes(x = year+2000, y= median, color = variable)) + 
  xlab("Year") + 
  ylab("Alpha")+
  theme_cowplot()

p2 <- alpha %>%
  ggplot() + 
  geom_line(aes(x = year+2000, y= median_transformed, color = variable)) + 
  xlab("Year") + 
  ylab(expression(paste(g^{-1}, "(mu+alpha)")))+
  theme_cowplot()

alpha_plot <- plot_grid(
  plot_grid(
    p1 + theme(legend.position="none"),
    p2 + theme(legend.position="none"),
    labels = "AUTO"
  ),
  get_legend(p1 + theme(legend.position="bottom")),
  ncol =1,
  rel_heights = c(5,1)
)

ggsave("alpha_plot.png", alpha_plot)



sigma_alpha <- phi_alpha_summary %>% filter(str_detect(variable, "sigma_alpha_pi")) %>% pull(median) %>% as.numeric()
sigma_alpha
sigma_delta <- phi_alpha_summary %>% filter(str_detect(variable, "sigma_delta_pi")) %>% pull(median)%>% as.numeric()
sigma_delta
sigma_lambda <- phi_alpha_summary %>% filter(str_detect(variable, "sigma_convolved_lambda")) %>% pull(median) %>% as.numeric()
sigma_lambda
sigma_pi <- phi_alpha_summary %>% filter(str_detect(variable, "sigma_convolved_pi")) %>% pull(median) %>% as.numeric()
sigma_pi
rho <- phi_alpha_summary %>% filter(str_detect(variable, "rho")) %>% pull(median)%>% as.numeric()
rho


# Table 1: Convergence Diagnostics
bind_rows(
  pi_summary %>%
    summarize(param = "pi",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  lambda_summary %>%
    summarize(param = "lambda",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "phi_pi"))%>%
    summarize(param = "phi_pi",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "theta_pi"))%>%
    summarize(param = "theta_pi",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "theta_lambda"))%>%
    summarize(param = "theta_lambda",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "alpha_pi"))%>%
    summarize(param = "alpha_pi",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "alpha_lambda"))%>%
    summarize(param = "alpha_lambda",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "delta_pi"))%>%
    summarize(param = "delta_pi",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ), 
  phi_alpha_summary %>%
    filter(str_detect(variable, "delta_lambda"))%>%
    summarize(param = "delta_lambda",
              rhat_q5 = quantile(rhat, 0.05),
              rhat_q95 = quantile(rhat, 0.95),
              ess_q5 = quantile(ess, 0.05),
              ess_q95 = quantile(ess, 0.95)
              ) 
) %>% as.data.frame() %>%
  stargazer::stargazer(summary = FALSE, digits = 2)

# Overall spatial pattern
overall_shp <- br_shp %>%
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "phi_pi")) %>%
      mutate(phi_pi = as.numeric(median))  %>%
      select(phi_pi)
  ) %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "theta_pi"))  %>%
      mutate(theta_pi = as.numeric(median)) %>%
      select(theta_pi)
  ) %>%
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "theta_lambda"))  %>%
      mutate(theta_lambda = as.numeric(median)) %>%
      select(theta_lambda)
  ) %>%
  mutate(across(starts_with("pop_"), as.numeric)) %>%
  mutate(pop = rowMeans(st_drop_geometry(select(., starts_with("pop_")))),
         pi= sigmoid(mu_pi + sigma_pi*(sqrt(rho)*phi_pi + sqrt(1-rho)*theta_pi)),
         lambda = exp(mu_lambda + sigma_lambda * theta_lambda),
         lambda_pop = as.numeric(pop) * exp(mu_lambda + sigma_lambda * theta_lambda),
         E_y= (1-pi) * lambda_pop)

# Overall spatial pattern
overall_plot <- tm_shape(overall_shp) + 
  tm_fill(col = "E_y", style="cont", title="Ey")  +
  tm_shape(ufs) +
  tm_borders(col="black") + 
  tm_text("abbrev_state", size = 0.5, col = "black") +
  tm_layout(main.title= 'Expected Cases: Overall Spatial Pattern', 
            main.title.size = 1)

tmap_save(overall_plot, filename = "overall_spatial.png", width=7, height = 5)

# Bernoulli process only
bern_plot <- tmap_arrange(
  tm_shape(overall_shp) +
    tm_fill(col = "pi", style = "cont", title = expression(pi)) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'A', 
              main.title.size = 1),
  
  tm_shape(overall_shp) +
    tm_fill(col = "phi_pi", style = "cont", title = expression(phi[pi])) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'B', 
              main.title.size = 1),
  
  tm_shape(overall_shp) +
    tm_fill(col = "theta_pi", style = "cont", title = expression(theta[pi])) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'C', 
              main.title.size = 1),
  nrow = 1
  
)

tmap_save(bern_plot, filename = "bern_plot.png", width = 7, height = 3)

# Lambda plot
lambda_plot <- tmap_arrange(
  tm_shape(overall_shp) +
    tm_fill(col = "lambda_pop", style = "cont", title = expression(Pop %*% lambda)) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'A', 
              main.title.size = 1),
  tm_shape(overall_shp) +
    tm_fill(col = "lambda", style = "cont", title = expression(lambda)) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'B', 
              main.title.size = 1),
  
  tm_shape(overall_shp) +
    tm_fill(col = "theta_lambda", style = "cont", title = expression(theta[lambda])) +
    tm_shape(ufs) +
    tm_borders(col="black") + 
    tm_text("abbrev_state", size = 0.5, col = "black") +
    tm_layout(main.title= 'C', 
              main.title.size = 1),
  nrow = 1
  
)

## Calculate Moran's I on theta
neighbours <- spdep::poly2nb(overall_shp)
listw <- spdep::nb2listw(neighbours, zero.policy = TRUE)
globalMoran <- spdep::moran.test(overall_shp$theta_lambda, listw, zero.policy = TRUE)
globalMoran

globalMoran <- spdep::moran.test(overall_shp$lambda, listw, zero.policy = TRUE)
globalMoran

tmap_save(lambda_plot, filename = "lambda_plot.png", width = 7, height = 3)
  




## From the climate model: 
## First fig: the estimated Betas
## Second fig: other parameter differences from main model

pi_summary_clim <- covariate_mcmc[["pi_summary"]] %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,",".")) %>%
  separate(variable, sep = "\\.", into = c("variable","year", "muni", NA))

lambda_summary_clim <- covariate_mcmc[["lambda_summary"]] %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,",".")) %>%
  separate(variable, sep = "\\.", into = c("variable","year", "muni", NA))

phi_alpha_summary_clim <- covariate_mcmc[["phi_alpha_summary"]] %>%
  mutate(across(c(median, rhat, ess), ~as.numeric(.))) 


E_y_clim <- pi_summary_clim  %>%
  mutate(variable = str_replace(variable, "pi", "E_y")) %>%
  select(variable, year, muni, median) %>% 
  mutate(year = as.numeric(year) + 2000)
E_y_clim$median  <- (1-pi_summary_clim$median)*lambda_summary_clim$median
E_y_clim
E_y_clim <- E_y_clim %>% left_join(muni_crosswalk, by = c("muni", "year"))
E_y


## Calculate RMSE
clim_rmse <- left_join(
  E_y_clim,
  br_shp %>%
    st_drop_geometry() %>% 
    select(muni_code, starts_with("chagas")) %>%
    pivot_longer(-muni_code, values_to = "chagas") %>% 
    separate(name, into = c(NA, "year"), convert=TRUE) 
)

clim_rmse %>% 
  # group_by(year) %>%
  summarize(rmse = sqrt(mean((median - chagas)^2)))


## Need to get distribution of Betas from the draws file
betas <- readRDS("covariate_phi_alpha_draws.RDS")
mus <- posterior::as_draws(map(betas, ~.[str_detect(names(.), "mu")]))
betas <- posterior::as_draws(map(betas, ~.[str_detect(names(.), "beta")]))
betas_summary <- summary(betas)
mus_summary <- summary(mus)
bioclimate_pc <- readRDS("bioclimate_pc.RDS")

## Plot of betas
betas_plot <- bayesplot::mcmc_areas(betas, regex_pars = c("^beta_pi[0-9]*", "^beta_lambda[0-9]*")) + 
  cowplot::theme_cowplot() + 
  geom_vline(aes(xintercept = 0), color = "red")
ggsave("betas_plot.png", betas_plot, height = 5, width = 7)


beta_original <- bind_cols(
  `pi (probability)` = sigmoid(mus_summary[[1,3]] + (as.matrix(bioclimate_pc$rotation[,1:6]) %*% as.vector(betas_summary %>% filter(str_detect(variable, "^beta_pi\\[*")) %>% pull(mean)))[,1]),
  `lambda (per million)` = 1000000*exp(mus_summary[[2,3]] + (as.matrix(bioclimate_pc$rotation[,1:6]) %*% as.vector(betas_summary %>% filter(str_detect(variable, "^beta_lambda\\[*")) %>% pull(mean)))[,1]),
  name = rownames(bioclimate_pc$rotation) 
)  


beta_original_plot <- ggplot(beta_original %>% pivot_longer(-name, names_to = "var")) + geom_tile(aes(x = var, y = name, fill = value)) + 
  viridis::scale_fill_viridis() +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("Process") +
  ggtitle("Transformed Coefficients for\n19 Bioclimatic Variables")

ggsave("beta_original_plot.png", beta_original_plot, height = 7, width = 8)


## Parameter differences from main model
model_diffs <- left_join(
  phi_alpha_summary, phi_alpha_summary_clim, suffix = c(".main", ".clim"), 
  by = c("variable"))

model_diffs <- model_diffs %>% mutate(diff = median.main - median.clim,
                                      variable = str_replace(variable, "\\[[0-9+,]*\\]", ""))

model_diffs_plot <- cowplot::plot_grid(
model_diffs %>%
  group_by(variable) %>% 
  filter(n() == 1) %>%
  ggplot() + geom_col(aes(variable, diff)) + 
  cowplot::theme_cowplot() + 
  xlab("")+
  ylab("") + 
  coord_flip() + 
  ggtitle("Parameter Differences between main and climate model"),
model_diffs %>%
  group_by(variable) %>% 
  filter(n() > 1) %>%
  ungroup() %>%
  ggplot() + geom_histogram(aes(diff)) + facet_wrap(~variable, scales = "free_y", nrow = 1) + 
  xlab("") + 
  ylab("") +
  cowplot::theme_cowplot(),

ncol = 1 

)

ggsave("model_diffs_plot.png", model_diffs_plot)


## Look at phi in the climate model
phi_clim_plot <- bind_cols(br_shp,
  phi_clim = model_diffs %>%  
      filter(str_detect(variable, "phi_pi")) %>% pull("diff")) %>%
  tm_shape() + 
  tm_fill(col = "phi_clim", style = "cont", title = expression(phi[pi]))+
  tm_shape(ufs) +
  tm_borders(col="black") + 
  tm_text("abbrev_state", size = 0.5, col = "black") +
  tm_layout(main.title= 'B', 
            main.title.size = 1)
  

tmap_save(phi_clim_plot, filename = "phi_clim_plot.png")

stop("The rest are old figs, not used")

## The rest are old figs, not used




# Plot phi and theta
phi_shp <- br_shp %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "phi_pi")) %>%
      mutate(median_transformed = sigmoid(mu_pi + median))
  )
theta_shp <- br_shp %>% 
  bind_cols(
    phi_alpha_summary %>% filter(str_detect(variable, "theta_pi"))  %>%
      mutate(median_transformed =  sigmoid(mu_pi + median))
  )

br_phi_theta <- tmap_arrange(
  tm_shape(phi_shp) + 
    tm_fill(col = "median", style="cont", title="Phi")  +
    tm_shape(ufs) +
    tm_borders() + 
    tm_text("abbrev_state", size = 0.5) +
    tm_layout(main.title= 'Pr(pi=0): Clustering term', 
              main.title.size = 1) , 
  tm_shape(theta_shp) + 
    tm_fill(col = "median", style="cont", title="Theta") +
    tm_shape(ufs) +
    tm_borders() + 
    tm_text("abbrev_state", size = 0.5) +
    tm_layout(main.title= 'Pr(pi=0):  Heterogeneity Term',
              main.title.size = 1),
  nrow=1
)

# tmap_mode("plot")
# br_phi_theta

tmap_save(br_phi_theta, filename = "fig2.png", width=7, height = 5)



# Group E_y by state
E_y_summary <- E_y_summary %>% 
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,", ".")) %>%
  separate(variable, sep="\\.", into = c("variable", "year", "muni"))

E_y_shp <- E_y_summary %>%
  select(year, muni, median) %>%
  pivot_wider(names_from = year, values_from = median, names_prefix = "year_")

E_y_shp <- bind_cols(br_shp, E_y_shp)

E_y_shp <- E_y_shp %>%
  rowwise() %>%
  mutate(populacao = as.numeric(populacao)) %>%
  mutate(across(contains("year_"), ~./populacao)) %>%
  mutate(sum_E = sum(c_across(starts_with("year"))))

# Calculate state sum
E_y_shp %>%
  st_drop_geometry() %>%
  group_by(uf) %>% 
  summarize(total_E_y = sum(sum_E)) %>%
  arrange(-total_E_y)

# Possible that theta and phi are not identified, or that theta is dominating
# over phi. Possible that BYM2-convolution can help...
E_y_shp %>%
  st_drop_geometry() %>%
  filter(uf == "MG") %>%
  pivot_longer(contains("year"), names_to = "year") %>%
  mutate(year = as.numeric(str_replace(year, "year_", ""))) %>%
  ggplot() + 
    geom_line(aes(x = year, y = value, group=nome))


# Alpha over time 
alpha_plot <- phi_alpha_summary %>%
  filter(str_detect(variable, "alpha")) %>%
  filter(!str_detect(variable, "sigma_alpha_pi")) %>%
  mutate(variable = str_replace_all(variable, "\\[|\\]|\\,", ".")) %>%
  separate(variable, sep="\\.", into = c("variable", "year", NA)) %>%
  mutate(year = as.numeric(year)) %>%
  ggplot() + 
    geom_line(aes(2001:2019, sigmoid(mu + sigma_alpha*median))) + 
    xlab("Year")  + 
    ylab("Probability") + 
    theme_cowplot()

alpha_plot

ggsave("alpha_plot.png", alpha_plot, width=7, height=3)
# Relationships with naive IR


# E_y
E_y_summary

hist(E_y_summary$median)





ggplot() + geom_point(aes(
  phi_shp$ess_bulk,
  theta_shp$ess_bulk
))

# Check pairs plot
ggplot() + geom_point(
  aes(
    extract_variable(phi_alpha, "phi_pi[1000]"),
    extract_variable(phi_alpha, "theta_pi[1000]")
  )
)

# Check convergence
ggplot() + geom_histogram(aes(E_y_summary$rhat)) + geom_vline(aes(xintercept = 1.01))
ggplot() + geom_histogram(aes(phi_alpha_summary$rhat)) + geom_vline(aes(xintercept = 1.01))

mean(E_y$ess_bulk)
mean(phi_alpha_summary$ess_bulk)

# Some variables clearly didn't converge:
E_y_summary %>%
  arrange(-rhat) %>% 
  slice_head(prop=0.01)

E_y_summary %>%
  arrange(ess_bulk) %>% 
  slice_head(prop=0.01)

chagas_summary %>% summarize(rhat_below_1.01 = mean(rhat < 1.01),
                             rhat_below_1.1 = mean(rhat < 1.1),
                             mean_ESS = mean(ess_bulk)) 


# Pick a few municipalities and plot their parameters over time
chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median))

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,190\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median)) 

chagas_summary %>% 
  filter(str_detect(variable, "^E_y\\[[0-9]*,193\\]")) %>%
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]") %>%
  ggplot() + geom_point(aes(x = t, y = median))

# Test mapping
IR <- chagas_summary %>%
  filter(str_detect(variable, "^E_y\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")
IR <- IR %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

# count <- as.data.frame(t(chagas_offset_count[1:5, ])) 
# colnames(count) = c("MLE_y1","MLE_y2", "MLE_y3", "MLE_y4", "MLE_y5")
# 
# IR <- bind_cols(IR, count) 
IR <- bind_cols(br_shp, IR)

IR <- IR %>% 
  mutate(populacao = log(as.numeric(populacao))) %>%
  mutate(across(contains("year"), ~./populacao))


tm_shape(IR) + tm_polygons(col = "year4", style = "cont", id = "nome")





pi <- chagas_summary %>%
  filter(str_detect(variable, "^pi\\[")) %>% 
  # split variable indices into columns
  separate(variable, into=c("variable", "t", "muni", NA), sep = "\\[|,|\\]")

pi <- pi %>%
  select(variable, t, muni, median) %>%
  pivot_wider(id_cols = c("variable","muni"), names_from = "t", names_prefix = "year", values_from = "median")

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