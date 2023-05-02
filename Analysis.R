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
load("./from_ayesha/cleaned_data/pop_all.Rdata")

to_drop <- br_shp %>% slice(1524, 3498, 5564) %>% pull(muni_code)
br_shp <- br_shp %>%filter(!muni_code %in% to_drop) %>% arrange(muni_code)

# st_write(br_shp, dsn = "br_shp.shp")


# Get map of UFs
ufs <- geobr::read_state()

logit <- function(x) {
  log(x/(1-x))
}

sigmoid <- function(x) {
  1/(1+exp(x))
}
# chagas_sample <- readRDS("mcmc_out/chagas_offset.RDS")

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


load <- FALSE
 
fname <- "mcmc_out/knorr_held_convolved_both_parts-202304281615-X-36f921.csv"
f1 <- str_replace(fname, "X", "1")
f2 <- str_replace(fname, "X", "2")
f3 <- str_replace(fname, "X", "3")
f4 <- str_replace(fname, "X", "4")


# Predicted counts

if (load) {
  load("pi_summary.Rdata")
  load("pi_draws.Rdata")
} else {
  pi1 <- read_cmdstan_csv(f1,variables = "pi", format = "draws_list")
  pi2 <- read_cmdstan_csv(f2,variables = "pi", format = "draws_list")
  pi3 <- read_cmdstan_csv(f3, variables = "pi", format = "draws_list")
  pi4 <- read_cmdstan_csv(f4, variables = "pi", format = "draws_list")
  
  pi <- bind_draws(x = pi1$post_warmup_draws, pi2$post_warmup_draws, 
                    pi3$post_warmup_draws, pi4$post_warmup_draws, along = "chain")
  rm(pi1, pi2, pi3, pi4)
  gc()
  
  # Since the datset is so big need to build the summary table by hand
  pi_names <- variables(pi)
  
  start.time <- Sys.time()
  pi_summary <- foreach(x = pi_names) %do% {
    d <- cbind(pi[[1]][[x]], pi[[2]][[x]], pi[[3]][[x]], pi[[4]][[x]])
    c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
  }
  pi_summary <- bind_rows(pi_summary)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken ## About 23 minutes
  
  save(pi_summary, file="pi_summary.Rdata")
  save(pi, file="pi_draws.Rdata", compress=TRUE)
  rm(pi)
}

if (load) {
  load("lambda_summary.Rdata")
  load("lambda_draws.Rdata")
} else {
  lambda1 <- read_cmdstan_csv(f1,
    variables = "lambda", format = "draws_list")
  lambda2 <- read_cmdstan_csv(f2,
    variables = "lambda", format = "draws_list")
  lambda3 <- read_cmdstan_csv(f3,
    variables = "lambda", format = "draws_list")
  lambda4 <- read_cmdstan_csv(f4,
    variables = "lambda", format = "draws_list")
  
  lambda <- bind_draws(x = lambda1$post_warmup_draws, lambda2$post_warmup_draws, 
                   lambda3$post_warmup_draws, lambda4$post_warmup_draws, along = "chain")
  rm(lambda1, lambda2, lambda3, lambda4)
  gc()
  
  # Since the datset is so big need to build the summary table by hand
  lambda_names <- variables(lambda)
  
  start.time <- Sys.time()
  lambda_summary <- foreach(x = lambda_names) %do% {
    d <- cbind(lambda[[1]][[x]], lambda[[2]][[x]], lambda[[3]][[x]], lambda[[4]][[x]])
    c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
  }
  lambda_summary <- bind_rows(lambda_summary)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken ## About 23 minutes
  
  save(lambda_summary, file="lambda_summary.Rdata")
  save(lambda, file="lambda_draws.Rdata", compress=TRUE)
}

# phi, alpha, theta, mu
if (load) {
  load("phi_alpha_summary.Rdata")
  load("phi_alpha_draws.Rdata")
} else {
  
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
    "sigma_delta_pi"
  )
  # Repeat for Phi, alpha, and theta, which are much samller 
  phi_alpha1 <- read_cmdstan_csv(f1,
    variables = other_vars, format = "draws_list")
  phi_alpha2 <- read_cmdstan_csv(f2,
    variables =  other_vars, format = "draws_list")
  phi_alpha3 <- read_cmdstan_csv(f3,
    variables =  other_vars, format = "draws_list")
  phi_alpha4 <- read_cmdstan_csv(f4,
    variables =  other_vars, format = "draws_list")
  
  phi_alpha <- bind_draws(x = phi_alpha1$post_warmup_draws, phi_alpha2$post_warmup_draws, 
                    phi_alpha3$post_warmup_draws, phi_alpha4$post_warmup_draws, along = "chain")
  
  phi_alpha_names <- variables(phi_alpha)
  phi_alpha_summary <- foreach(x = phi_alpha_names) %do% { 
    d <- cbind(phi_alpha[[1]][[x]], phi_alpha[[2]][[x]], phi_alpha[[3]][[x]], phi_alpha[[4]][[x]])
    c(variable = x, median = median(d), rhat = rhat(d), ess= ess_bulk(d))
  }
  phi_alpha_summary <- bind_rows(phi_alpha_summary)
  
  save(phi_alpha_summary, file= "phi_alpha_summary.Rdata")
  save(phi_alpha, file= "phi_alpha_draws.Rdata")

}


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
  select(nome, uf, muni_code, populacao) %>%
  mutate(populacao = as.numeric(populacao)) %>%
  rownames_to_column(var = "muni")


# Generate E_y from pi and Lambda
E_y  <- pi_summary  %>%
  mutate(variable = str_replace(variable, "pi", "E_y")) %>%
  select(variable, year, muni, median)
E_y$median  <- (1-pi_summary$median)*lambda_summary$median
E_y <- E_y %>% left_join(muni_crosswalk, by = c("muni"))
E_y

##Calculate incidence rate as sum(E_y)/sum(pop)

total_IR <- E_y %>% 
  summarize(IR = 100000*sum(median)/sum(populacao))
total_IR 

municipality_IR <- E_y %>%
  group_by(muni_code, uf, nome) %>%
  summarize(IR = sum(median)/sum(populacao))

municipality_IR

UF_IR <- E_y %>%
  group_by(uf) %>%
  summarize(IR = sum(median)/sum(populacao))

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
  mutate(year = as.numeric(year) + 1999)  

# IR
PA_AP %>% 
  group_by(nome) %>% 
  summarize(IR = 100000*sum(median)/sum(populacao),
            pop = mean(populacao)) %>%
  arrange(-IR) %>% View()

PA_AP %>% 
  group_by(nome) %>% 
  summarize(IR = sum(median),
            pop = mean(populacao))%>%
  arrange(-IR) %>% View()


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
  geom_line(aes(x = year+1999, y= median, color = variable)) + 
  xlab("Year") + 
  ylab("Alpha")+
  theme_cowplot()

p2 <- alpha %>%
  ggplot() + 
  geom_line(aes(x = year+1999, y= median_transformed, color = variable)) + 
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
  mutate(pi= sigmoid(mu_pi + sigma_pi*(sqrt(rho)*phi_pi + sqrt(1-rho)*theta_pi)),
         lambda = exp(mu_lambda + sigma_lambda * theta_lambda),
         lambda_pop = as.numeric(populacao) * exp(mu_lambda + sigma_lambda * theta_lambda),
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
    geom_line(aes(2000:2018, sigmoid(mu + sigma_alpha*median))) + 
    xlab("Year")  + 
    ylab("Probability") + 
    theme_cowplot()

alpha_plot

ggsave("alpha_plot.png", alpha_plot, width=7, height=5)
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