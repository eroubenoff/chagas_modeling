# TSIR for Brazil at the regional level
library(lfe)
library(dplyr)
library(lubridate)
library(MASS)

rm(list = ls())

########################
# get region codes
load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata") 
########################

# get dengue time series data by region
load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic_HOSPITAL.Rdata")
dat <- tmp7 %>% dplyr::select(code, den, month, year) %>%
  left_join(., region_master, by = "code") %>%
  group_by(region, year, month) %>%
  summarize(den = sum(den, na.rm = T)) %>%
  ungroup()
save(dat, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_region_tsir.Rdata")

########################
# get population time series by region (TCU estimate)
text <- read.csv("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/popTCUEst_1992-2017.csv", sep = ";", header = FALSE, fileEncoding="latin1")

# get column names
cols <- as.vector(apply(text[4,],1, function(x) as.character(x)))
cols[1] <- "name"

# get table and add column names
pop <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
colnames(pop) <- cols
pop$code <- as.numeric(substr(pop[,1],1,6))
pop[,2:ncol(pop)] <- sapply(pop[,2:ncol(pop)], as.character)
pop[,2:ncol(pop)] <- sapply(pop[,2:ncol(pop)], as.numeric)


# reshape to long format
keycol <- "year"
valuecol <- "ps"
gathercols <- colnames(pop)[which(colnames(pop) != "name" & colnames(pop) != "code")]
pop.all <- gather_(pop, keycol, valuecol, gathercols, na.rm = FALSE)
pop.all$code <- as.numeric(pop.all$code)
pop.all$ps <- as.numeric(as.character(pop.all$ps))

ps <- pop.all %>% left_join(.,region_master, by = "code") %>%
  group_by(region, year) %>%
  summarize(ps = sum(ps, na.rm = T)) %>%
  ungroup()
ggplot(ps) + geom_line(aes(x = year, y = ps, group = region))
save(ps, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/pop_region_tsir.Rdata")
########################
# get births time series by region

#1998-2001
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/Births")
filenames <- list.files(pattern = "Monthly")
full.births <- NA
for (i in 1:length(filenames)) {
  text <- read.csv(filenames[i], sep = ";", header = FALSE, fileEncoding="latin1")
  # get column names
  cols <- as.vector(apply(text[4,],1, function(x) as.character(x)))
  cols[1] <- "name"
  
  # get table and add column names
  tmp <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
  colnames(tmp) <- cols
  tmp$code <- as.numeric(substr(tmp[,1],1,6))
  tmp[,2:ncol(tmp)] <- sapply(tmp[,2:ncol(tmp)], as.character)
  tmp[,2:ncol(tmp)] <- sapply(tmp[,2:ncol(tmp)], as.numeric)
  
  # reshape to long format
  keycol <- "month"
  valuecol <- "births"
  gathercols <- colnames(tmp)[which(colnames(tmp) != "name" & colnames(tmp) != "code" & colnames(tmp) != "Total" & colnames(tmp) != "Ignorado")]
  tmp.all <- gather_(tmp, keycol, valuecol, gathercols, na.rm = FALSE) %>%
    dplyr::select(code, month, births) %>%
    mutate(code = as.numeric(code), births = as.numeric(as.character(births)))
  tmp.all$births[is.na(tmp.all$births)] <- 0

  tmp2 <- tmp.all %>% left_join(.,region_master, by = "code") %>%
    group_by(region, month) %>%
    summarize(births = sum(births, na.rm = T)) %>%
    ungroup()
  tmp2$year <- as.numeric(substr(filenames[i],8,11))
  
  full.births <- rbind(full.births, tmp2)
  
}

full.births <- full.births %>% filter(!is.na(full.births$region)) %>%
  rename(month1 = month) 
full.births <- full.births %>% mutate(month = case_when(month1 == "Janeiro" ~ 1,
                                                        month1 == "Fevereiro" ~ 2,
                                                        month1 == "Março" ~ 3,
                                                        month1 == "Abril" ~ 4,
                                                        month1 == "Maio" ~ 5,
                                                        month1 == "Junho" ~ 6,
                                                        month1 == "Julho" ~ 7,
                                                        month1 == "Agosto" ~ 8,
                                                        month1 == "Setembro" ~ 9,
                                                        month1 == "Outubro" ~ 10,
                                                        month1 == "Novembro" ~ 11,
                                                        month1 == "Dezembro" ~ 12))





# add in 2003-2016
text <- read.csv("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/Births/births20032016.csv", header = FALSE, fileEncoding="latin1")
tab <- text[8:5577, 2:ncol(text)] %>% mutate(code = as.numeric(as.character(V2))) %>%
  dplyr::select(code, everything()) %>% dplyr::select(-c(V2, V3, V4, V5)) %>%
  mutate(code = as.numeric(substr(as.character(code),1,6)))
  

years = 2003:2016
colnum = 2:13
full.births2 <- NA
for (y in 1:length(years)) {
  tmp <- tab[,c(1,colnum)]
  head(tmp)
  names(tmp)[2:13] <- 1:12
  keycol <- "month"
  valuecol <- "births"
  gathercols <- colnames(tmp)[which(colnames(tmp) != "code")]
  tmp.all <- gather_(tmp, keycol, valuecol, gathercols, na.rm = FALSE) %>%
    dplyr::select(code, month, births) %>%
    mutate(code = as.numeric(code), births = as.numeric(as.character(births)))
  tmp.all$births[is.na(tmp.all$births)] <- 0
  
  tmp2 <- tmp.all %>% left_join(.,region_master, by = "code") %>%
    group_by(region, month) %>%
    summarize(births = sum(births, na.rm = T)) %>%
    ungroup()
  tmp2$year <-years[y]
  
  full.births2 <- rbind(full.births2, tmp2)
  colnum = colnum + 12
  
  
}

full.births2 <- full.births2 %>% filter(!is.na(full.births2$region))
all.births <- full.births %>% dplyr::select(region, month, births, year) %>%
  rbind(.,full.births2)
head(all.births)



save(all.births, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/births_region_tsir.Rdata")


########################

### Simulation function
SimTsir <- function(ppsize,
                    startS,
                    startI,
                    beta_ts = beta_ts,
                    alpha=0.97,
                    births,
                    m=0.1){
  
  Ival <- Sval <- extinct <- rep(NA,length(births));
  Ival[1] <- startI;
  Sval[1] <- startS
  
  # calculate beta for each time step
  beta = beta_ts
  
  for (t in 1:(length(beta)-1)){
    
    
    #1. use neg binomial Xia style
    muval <- (beta[t]*(Ival[t]^alpha)*Sval[t])/ppsize[t]
    Ival[t+1] <- rnegbin(1,muval,pmax(Ival[t],1))
    
    #2. Alternative + simpler stoch bit! hash out previous line and unhash this one to explore...
    #Ival[t+1] <- rbinom(1,Sval[t],1-exp(-(beta[season[t]]*(Ival[t]^alpha))/ppsize[t]))
    #hash out or not depending on how you feel about this - definitely helps get over transients at start....
    if (Ival[t+1]==0) {extinct[t+1]=1; Ival[t+1]=1}
    
    #immigrants
    Ival[t+1] <- Ival[t+1] + rpois(1,m)
    
    #update susceptibles
    Sval[t+1] <- max(Sval[t]-Ival[t+1]+births[t],0) #can't have negative susceptibles!
    print(t)
  }
  
  return(list(Ival=Ival,Sval=Sval,extinct=extinct))
}




### Load data, interpolate, then do TSIR

rm(list = ls())
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/births_region_tsir.Rdata")
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/pop_region_tsir.Rdata")
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/dengue_region_tsir.Rdata")

for (r in 1:length(unique(dat$region))){
  
  time = seq(2003, 2017, by = 1/24)[1:((length(seq(2003, 2017, by = 1/24)) - 1))]
  
  # cases
  cases.month <- dat %>% filter(region == r & as.numeric(as.character(year)) >= 2003 & as.numeric(as.character(year)) <= 2016) %>%
    arrange(year, month) %>%
    mutate(time = seq(2003, 2017, by = 1/12)[1:nrow(.)])
  
  cases <- round(approx(cases.month$time,cases.month$den, xout = time, 
                        rule = 2, method = "linear", ties = mean)$y * 12/24)
  
  #births
  births.month <- all.births %>% filter(region == r & as.numeric(as.character(year)) >= 2003 & as.numeric(as.character(year)) <= 2016) %>%
    arrange(year, month) %>%
    mutate(time = seq(2003, 2017, by = 1/12)[1:nrow(.)])

  births <- round(approx(births.month$time,births.month$births, xout = time, 
                        rule = 2, method = "linear", ties = mean)$y * 12/24)
  
  #population
  pop.year <- ps %>% filter(region == r & as.numeric(as.character(year)) >= 2003 & as.numeric(as.character(year)) <= 2016) %>%
    arrange(year) %>%
    mutate(time = as.numeric(as.character(year)))
  
  pop <- round(approx(pop.year$time,pop.year$ps, xout = time, 
                      rule = 2, method = "linear", ties = mean)$y)
  
  # now try TSIR
  lengthdata <- length(time)
  
  B <- births
  C <- cases
  
  X <- cumsum(C) # cumulative cases
  Y <- cumsum(B) # cumulative births
  
  cum.reg<-smooth.spline(Y,X,df=4) #Regressing X on Y
  ur <- predict(cum.reg,Y,deriv=1)$y #ur = 1/rho
  Z <- (predict(cum.reg)$y - X) * (1/ur)
  
  m.rho <- mean(ur)
  
  # Ic are actual cases: reported cases multiplied by rho
  Ic= C/ur # Ic = r
  plot(Ic, type = "l")
  
  
  ## Keep only positive values
  keep.index <- which(Ic[1:(lengthdata - 1)] > 0 & Ic[2:lengthdata] > 0)
  length.pos <- length(keep.index)
  
  lInew <- log(Ic[2:lengthdata])[keep.index]
  Inew <- exp(lInew)
  lIold <- log(Ic[1:(lengthdata -1)]) [keep.index]
  Zold = Z[1:(lengthdata - 1)] [keep.index]
  
  year <- substr(time, 1,4)
  table(year)
  
  seas <- rep(1:24, length(unique(year)))
  seas <- seas[2:lengthdata][keep.index]
  
  ps = pop[1:(lengthdata - 1)] [keep.index]
  
  form = round(Inew) ~ -1 + as.factor(seas) + offset(0.97*lIold) + offset(lSold - log(ps))
  
  
  Smean = seq(0.02, 0.5, by = 0.01)
  llik = rep(NA, length(Smean))
  
  for (i in 1:length(Smean)) {
    lSold = log(Smean[i]*mean(ps) + Zold)
    
    glmfit = glm(formula = form, family = poisson(link = "log"))
    
    llik[i] = glmfit$deviance # -2*loglikelihood
  }
  
  Sbar <- Smean[which(llik == min(llik))] #we want to maximize the log likelihood
  
  plot(Smean, llik, type = "l")
  abline(v = Sbar)
  
  
  lSold = log(Sbar*mean(ps) + Zold)
  # best estimates
  lmfit = glm(formula = form, family = poisson(link = "log"))
  sum <- summary(lmfit)
  
  beta.tmp = exp(lmfit$coeff[1:length(unique(seas))])
  indexes <- as.numeric(substring(names(beta.tmp),16,nchar(names(beta.tmp))))
  beta = rep(NA, max(seas))
  beta[indexes] = beta.tmp  
  beta_ts = beta[seas]
  beta_ts[is.na(beta_ts)] <- mean(beta, na.rm = TRUE)
  
  plot(beta, type = "l")
  
  testSim <- SimTsir(ppsize = pop, startS = exp(lSold[1]), startI = exp(lIold[1]),
                     beta_ts = beta_ts, alpha = 0.97, births = B)
  plot(Ic, type = "l")
  lines(testSim$Ival, col = "red")
  
}



########################
