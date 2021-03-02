library(foreign)
library(dplyr)
library(stringr)
library(gdata)
library(pracma)
rm(list = ls())

#############################################################################
# get pop data and names and codes for all municipalities
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData")
text<- read.csv(list.files()[1], sep = ";", header = FALSE, fileEncoding="latin1")

# get column names
cols <- as.vector(apply(text[4,],1, function(x) as.character(x)))
cols[1] <- "name"

# get table and add column names
pop <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
colnames(pop) <- cols
pop$code <- as.numeric(substr(pop[,1],1,6))

save(pop, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")

#############################################################################
#get full list of municipality names and codes
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")
states <- pop[,c("name","code")]

#############################################################################
# clean dengue data
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = "csv")

for (i in 1:length(filenames)) {
  text <- read.csv(filenames[i], sep = ";", header = FALSE, fileEncoding="latin1")
  
  #get year
  year =  substr((as.character(text[3,1])), (regexpr(':', text[3,1])[1]+1), nchar(as.character(text[3,1])))
  
  # get column names
  cols <- as.character(unlist(text[4,]))
  cols[1] <- "name"
  cols <- paste(cols, "_d",sep = "")
  
  
  # get table and add column names
  tab <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
  colnames(tab) <- cols
  tab$code <- as.numeric(substr(tab[,1],1,6))
  
  #drop column on missing values
  if (colnames(tab)[2] == "Em branco/IGN_d" | colnames(tab)[2] == "Ign/Em Branco_d") tab <- tab[,-2]
  
  #turn state names to character  
  tab[,1] <- as.character(tab[,1])
  
  
  # convert case count to numeric, assign the zero values for months where there were no cases
  tab[,2:which(colnames(tab) == "Total_d")] <- sapply(tab[,2:which(colnames(tab) == "Total_d")],as.character)
  tab[which(tab == "-", arr.ind = TRUE)] <- "0"
  tab[,2:which(colnames(tab) == "Total_d")] <- sapply(tab[,2:which(colnames(tab) == "Total_d")],as.numeric)
  
  data <- left_join(states, tab, by = "code")
  data[is.na(data)] = 0
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",year,".Rdata", sep = "")
  save(data, file = filepath)
  print(i)
  
}
#############################################################################
# clean malaria data
#setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN")

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion")
filenames <- list.files(pattern = "csv")

for (i in 1:length(filenames)) {
  text <- read.csv(filenames[i], sep = ";", header = FALSE, fileEncoding="latin1")
  
  #get year
  year =  substr((as.character(text[3,1])), (regexpr(':', text[3,1])[1]+1), nchar(as.character(text[3,1])))
  
  # get column names
  cols <- as.character(unlist(text[4,]))
  cols[1] <- "name"
  cols <- paste(cols, "_m",sep = "")
  
  
  # get table and add column names
  tab <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
  colnames(tab) <- cols
  tab$code <- as.numeric(substr(tab[,1],1,6))
  
  #drop column on missing values
  if (colnames(tab)[2] == "Em branco/IGN_m" | colnames(tab)[2] == "Ign/Em Branco_m") tab <- tab[,-2]
  
  #turn state names to character  
  tab[,1] <- as.character(tab[,1])
  
  
  # convert case count to numeric, assign the zero values for months where there were no cases
  tab[,2:which(colnames(tab) == "Total_m")] <- sapply(tab[,2:which(colnames(tab) == "Total_m")],as.character)
  tab[which(tab == "-", arr.ind = TRUE)] <- "0"
  tab[,2:which(colnames(tab) == "Total_m")] <- sapply(tab[,2:which(colnames(tab) == "Total_m")],as.numeric)
  
  data <- left_join(states, tab, by = "code")
  data[is.na(data)] = 0
  
  
  filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Malaria-SINAN-infectionregion/MalariaSINAN_",year,".Rdata", sep = "")
  save(data, file = filepath)
  print(i)
  
}
#############################################################################
# get urban-rural pop percentage
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/Census2010/Brasil")
dat <- read.csv("tab1.csv", header = FALSE)
tmp <- dat[10:5580, c(1,4:6)] 
names(tmp) <- c("code_long", "pop2010", "Urban", "Rural")
tmp$code <- as.numeric(substr(tmp$code_long,1,6))
urban <- tmp[-c(which(is.na(tmp$code) == TRUE)),]

save(urban, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned.Rdata")


#urban pop 2000 (not %)
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/")
dat <- read.csv("urbanPop1991_2000.csv", header = FALSE)
tmp <- dat[2:5597, c(2:5)] 
names(tmp) <- c("code_long", "munic", "urban1991", "urban2000")
tmp$code <- as.numeric(substr(tmp$code_long,1,6))
urban <- tmp
save(urban, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_cleaned_1991-2000.Rdata")

## create urban-rural yearly %
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/")
dat.u <- read.csv("urbanPop_by_munic.csv", header = TRUE)
dat.r <- read.csv("ruralPop_by_munic.csv", header = TRUE)

dat.u <- dat.u[order(dat.u$Code), ]
dat.r <- dat.r[order(dat.r$Code), ]

perc1996 <- dat.u$X1996/(dat.u$X1996 + dat.r$X1996)
perc2000 <- dat.u$X2000/(dat.u$X2000 + dat.r$X2000)
perc2007 <- dat.u$X2007/(dat.u$X2007 + dat.r$X2007)
perc2010 <- dat.u$X2010/(dat.u$X2010 + dat.r$X2010)

dat <- data.frame(cbind(dat.u$Code, perc1996, perc2000, perc2007, perc2010))

interpolate <- function(x) {
  perc <- as.vector(unlist(x))
  if(is.na(perc[1])) {
    rep(NA, length = length(c(2000:2010)))
  } else {
    pracma::interp1(c(2000,2007,2010), perc, c(2000:2010), method = "linear")
    
  }
  
}


interped <- t(apply(dat[,c(3:5)],1, function(x) interpolate(x)))
urb <- data.frame(cbind(dat$V1, interped))
names(urb) <- c("code_long", paste0(2000:2010))
urb$code <- as.numeric(substr(urb$code_long,1,6))

save(urb, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_interp_2000-2010.Rdata")


# extrapolate
library(Hmisc)
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/")
dat.u <- read.csv("urbanPop_by_munic.csv", header = TRUE) %>% dplyr::select(Code, X1996, X2000, X2007, X2010)
dat.r <- read.csv("ruralPop_by_munic.csv", header = TRUE) %>% dplyr::select(Code, X1996, X2000, X2007, X2010)

dat.u <- dat.u[order(dat.u$Code), ]
dat.r <- dat.r[order(dat.r$Code), ]

perc1996 <- dat.u$X1996/(dat.u$X1996 + dat.r$X1996)
perc2000 <- dat.u$X2000/(dat.u$X2000 + dat.r$X2000)
perc2007 <- dat.u$X2007/(dat.u$X2007 + dat.r$X2007)
perc2010 <- dat.u$X2010/(dat.u$X2010 + dat.r$X2010)

dat <- dat.u %>% left_join(., dat.r, by = "Code")
names(dat)[1] <- "code_long"

dat2 <- dat %>% left_join(., region_master, by = "code_long") %>%
  group_by(region) %>%
  summarize(X1996.u = sum(X1996.x, na.rm = T),
            X1996.r = sum(X1996.y, na.rm = T),
            X2000.u = sum(X2000.x, na.rm = T),
            X2000.r = sum(X2000.y, na.rm = T),
            X2007.u = sum(X2007.x, na.rm = T),
            X2007.r = sum(X2007.y, na.rm = T),
            X2010.u = sum(X2010.x, na.rm = T),
            X2010.r = sum(X2010.y, na.rm = T)) %>%
  ungroup() %>% na.omit(.) %>%
  mutate(perc1996 = X1996.u/(X1996.u + X1996.r),
         perc2000 = X2000.u/(X2000.u + X2000.r),
         perc2007 = X2007.u/(X2007.u + X2007.r),
         perc2010 = X2010.u/(X2010.u + X2010.r))




interpolate <- function(x) {
  perc <- as.vector(unlist(x))
  #pracma::interp1(c(1996,2000,2007,2010), perc, c(1996:2017), method = "linear")
  Hmisc::approxExtrap(c(1996,2000,2007,2010), perc, c(1996:2017), method = "linear")$y
  
  
}



interped <- t(apply(dat2[,c(10:13)],1, function(x) interpolate(x)))
urb <- data.frame(cbind(dat2$region, interped))
names(urb) <- c("region", paste0(1996:2017))

save(urb, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/urban_interp_1996-2017.Rdata")





#############################################################################
# get land area
# data from https://ww2.ibge.gov.br/english/geociencias/cartografia/default_territ_area.shtm
# area is in km2
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData")
dat <- read.csv("area-by-munic.csv", header = TRUE)
tmp <- dat[1:5572, c(5,7)]
names(tmp) <- c("code_long", "area_km2")
tmp$code <- as.numeric(substr(tmp$code_long,1,6))
area <- tmp
save(area, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")

#############################################################################
# create file of municipalities by state and region

setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData")
dat <- read.csv("area-by-munic.csv", header = TRUE)
tmp <- dat[1:5572, c("CD_GCUF", "CD_GCMUN")]
names(tmp) <- c("stateID","code_long")
tmp$code <- as.numeric(substr(tmp$code_long,1,6))
tmp$stateID <- as.numeric(as.character(tmp$stateID))

#Need to manually code in region in case we want to aggregate up
#North Region = 1
## States: Acre, Amapá, Amazonas, Pará, Rondônia, Roraima, Tocantins
NR <- c(12, 16, 13,15, 11, 14, 17)

# Northeast Region = 2
## States: Alagoas, Bahia, Ceará, Maranhão, Paraíba, Pernambuco, Piauí, Rio Grande do Norte, Sergipe
NER <- c(27, 29, 23, 21, 25, 26, 22, 24, 28)

#Central West Region = 3
## States: Goiás, Mato Grosso, Mato Grosso do Sul, Distrito Federal (Federal District).
CWR <- c(52, 51, 50,53)

# Southeast Region = 4
## States: Espírito Santo, Minas Gerais, Rio de Janeiro, São Paulo
SER <- c(32, 31, 33, 35)

# South Region = 5
## States: Paraná, Rio Grande do Sul, Santa Catarina
SR <- c(41, 43, 42)


tmp$region <- ifelse(tmp$stateID %in% NR, 1,
                      ifelse(tmp$stateID %in% NER, 2,
                             ifelse(tmp$stateID %in% CWR, 3,
                                    ifelse(tmp$stateID %in% SER, 4,
                                           ifelse(tmp$stateID %in% SR, 5, NA)))))

region_master <- tmp

# add in state names
setwd("~/Google Drive/Pertussis")
text <- read.csv("ALL.csv", sep = ";", header = FALSE, fileEncoding="latin1")
tab_full <- text[5:(max(which(text == "Total", arr.ind = T)[,1])-1),]
cols <- as.character(unlist(text[4,]))
colnames(tab_full) <- cols
tab_full <- tab_full[,-2]
tab_full$stateID <- as.numeric(substr(tab_full[,1],1,2))
tab_all_states <- tab_full[,c("UF de notificação", "stateID")]
numchar <- nchar(as.character(tab_all_states[,1]))
for (x in 1:nrow(tab_all_states)) {
  tab_all_states$stateName[x] <- substr(as.character(tab_all_states[x,1]),4,numchar[x])
}
names.tmp <- tab_all_states[,c("stateName", "stateID")]

region_master <- left_join(region_master, names.tmp, by = "stateID")

save(region_master, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")

#############################################################################
# get GDP per capita
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/")
dat <- read.csv("GdpPerCapita2000-2012.csv", header = TRUE)
dat$code <- as.numeric(substr(dat[,1],1,6))
gdp <- dat[,-which(colnames(dat) == "Total")]
gdp <- gdp[,-which(colnames(dat) == "Year")]
save(gdp, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/gdp_cleaned.Rdata")

# agric value
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/")
dat <- read.csv("GDPAgricValue2000-2012.csv", header = TRUE)
dat$code <- as.numeric(substr(dat[,1],1,6))
agric <- dat[,-which(colnames(dat) == "Year")]
save(agric, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/agric_cleaned.Rdata")


#############################################################################
