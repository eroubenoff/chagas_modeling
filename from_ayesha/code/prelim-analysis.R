library(dplyr)
library(signal)
require(bbmle)

############### functions for estimating the attack rate
likelihood = function(S0, beta, I) {
  n = length(I)
  S = floor(S0 - cumsum(I[-n]))
  p = 1 - exp(-beta * (I[-n])/S0)
  L = -sum(dbinom(I[-1], S, p, log = TRUE))
  return(L)
}

sim.cb = function(I0,S0, beta) { 
  I=I0
  S = S0 
  i=1
  while (!any(I == 0)) { 
    i=i+1
    I[i] = rbinom(1, size = S[i - 1], prob = 1 - exp(-beta * I[i - 1]/S0))
    S[i] = S[i - 1] - I[i]
  }
  out = data.frame(S = S, I = I)
  return(out)
}

###############

############### make a map with before and after peak cutoffs for each state
bp.map <- ap.map <- list()
bp.map[[1]] <- c(16, 40, 64,  88, 112, 136, 160, 187, 208, 232, 256, 280)
ap.map[[1]] <- c(36,  60,  84, 108, 132, 156, 180, 206, 228, 252, 281)
bp.map[[2]] <- c(88, 114, 136, 160, 186, 208, 232, 256, 280)
ap.map[[2]] <- c(108, 124, 156, 180, 204, 228, 252, 276)
bp.map[[3]] <- c(18, 45, 64,  93, 163, 184, 236, 258, 280)
ap.map[[3]] <- c(36,  64,  84, 111, 180, 204, 252, 272)
bp.map[[4]] <- c(12, 58,  103, 160, 188, 214)
ap.map[[4]] <- c(30,  72, 112, 204, 200, 234)
bp.map[[5]] <- seq(17,288, by = 24)
ap.map[[5]] <- seq(36,288, by = 24)
###############

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/region_master_cleaned.Rdata")
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion")
filenames <- list.files(pattern = ".Rdata")
years <- unlist(lapply(filenames, function(x) substr(x,13,16)))
time <- seq(2001, 2013, 1/12)[1:144]
mtrx <- matrix(NA, 144, 27)
time.range <- seq(1,144, by = 12)

for (y in 1:length(years)){
  load(paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/Dengue-SINAN-infectionregion/DengueSINAN_",years[y],".Rdata", sep = ""))
  tmp <- data[,c(2,4:15)]
  tmp2 <- inner_join(region_master, tmp, by = "code")
  tmp2[is.na(tmp2)]<- 0
  df = tmp2[,c(1,6:17)] %>% group_by(stateID) %>% summarise_each(funs(sum))
  
  mtrx[time.range[y]:(time.range[y]+11),1:27] <- t(as.matrix(df[,2:13]))
}


#pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/dengue_timeseries_states2.pdf", width = 9, height = 9)
par(mfrow = c(5,6))
par(mar = c(4,3,3,2))
for (i in 1:27){
  plot(time, mtrx[,i], type = "l", main = unique(region_master$stateName)[i], xlab = "", ylab = "")
  
}
#dev.off()

i = 5
#interpolate cases
periodicity = 24
lengthdata <- length(time) * periodicity/12
time.bwk <- seq(2001,2013, 1/24)[1:lengthdata]
tstep <- as.numeric(as.factor(round((time.bwk %% 1), 3)))

cases <- as.numeric(as.vector(mtrx[,i]))
df <- as.data.frame(cbind(time, cases))

I <- round(interp1(df$time, df$cases, time.bwk, method = "linear", extrap = NA) * 12/periodicity)

#library(pracma)

plot(time.bwk, I, type = "l", main = unique(region_master$stateName)[i], xlab = "", ylab = "")
#x = findpeaks(I, nups = 2, ndowns = 2, minpeakheight = 200)
#x
#abline(v = time.bwk[x[,3]], lty =2, col = "red")
#abline(v = time.bwk[x[,4]], lty =2, col = "blue")
xseq <- seq(16,288, by = 24)
abline(v = time.bwk[xseq], lty = 2, col = "red")
yseq <- seq(36,288, by = 24)
abline(v = time.bwk[yseq], lty = 2, col = "blue")





ar <- ts.pk <- rep(NA, length = (length(xseq) -1 ))
ts <- ts.tstep <- list(NA)
for (pk in 1:(length(xseq)-1)) {
  I.tmp = I[xseq[pk]:yseq[pk]]
  ts.tmp = time.bwk[xseq[pk]:yseq[pk]]
  plot(ts.tmp, I.tmp, type = "l")
  fit = mle2(likelihood, start = list(S0 = 100000, beta = 3), method = "Nelder-Mead",
             data = list(I = I.tmp))
  summary(fit)
  plot(I.tmp, type = "n")
  for (x in 1:100) {
    sim = sim.cb(I0 = I.tmp[1],S0 = floor(coef(fit)["S0"]), beta = coef(fit)["beta"])
    lines(sim$I, col = grey(0.5))
  }
  points(I.tmp, type = "b", col = 2)
  ar[pk] = sum(I.tmp) / floor(coef(fit)["S0"])
  ts[[pk]] <- ts.tmp
  ts.tstep[[pk]] <- tstep[xseq[pk]:yseq[pk]]
  ts.pk[pk] <- ts.tmp[which(I.tmp == max(I.tmp, na.rm = TRUE))]
  
  
  
}



# get covariates
# population
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/pop_cleaned.Rdata")
pop_reg <- inner_join(region_master, pop, by = "code")

cols.num <- c(7:22)
pop_reg[cols.num] <- sapply(pop_reg[cols.num],as.numeric)
pop.df = pop_reg[,c(1,7:22)] %>% group_by(stateID) %>% summarise_each(funs(sum(., na.rm = TRUE)))
pop.tmp <- pop.df[pop.df$stateID == unique(region_master$stateID)[i],]


# area
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/PopData/area_cleaned.Rdata")
area$area_km2 <- as.numeric(as.character((area$area_km2)))
area.df <- inner_join(region_master, area, by = "code")
area.df2 = area.df[,c(1,7)] %>% group_by(stateID) %>% summarise_each(funs(sum(., na.rm = TRUE)))
area.tmp <- area.df2$area_km2[area.df2$stateID == unique(region_master$stateID)[i]]


ps <- dens <- rep(NA, length = (length(ar)))
for (t in 1:length(ar)) {
  pk <- round(ts.pk[t],0)
  ps[t] <- pop.tmp[,c(as.character(pk))]
  dens[t] <- pop.tmp[,c(as.character(pk))] / area.tmp 
}


# climate
clim <- list()
var = c("pre", "tmn", "tmp", "tmx")

for (t in 1:length(ar)) {
  min.yr <- floor(min(ts[[t]])); min.mth <- floor((ts.tstep[[t]][1])/2)
  max.yr <- floor(max(ts[[t]])); max.mth <- floor((ts.tstep[[t]][length(ts.tstep[[t]])])/2)

  clim.tmp <- rep(NA, length = length(var))
  for (v in 1:length(var)) {
    filepath <- paste("~/Google Drive/Dengue/Dengue-Malaria/Brazil/climdata/cleaned_",var[v],".Rdata", sep = "")
    load(filepath)
    data <- data.all
    
    min.col <- paste("Y", min.yr,"M",min.mth, sep = "")
    max.col <- paste("Y", max.yr,"M",max.mth, sep = "")
    data.tmp <- data[,c(which(colnames(data) == min.col):which(colnames(data) == max.col), ncol(data))]
    data.tmp[1:ncol(data.tmp)] <- sapply(data.tmp[1:ncol(data.tmp)],as.character)
    data.tmp[1:ncol(data.tmp)] <- sapply(data.tmp[1:ncol(data.tmp)],as.numeric)
    
    data.tmp2 <- inner_join(region_master, data.tmp, by = "code")
    
    data.tmp3 = data.tmp2[,c(1,6:ncol(data.tmp2))] %>% group_by(stateID) %>% summarise_each(funs(mean(., na.rm = TRUE)))
    data.tmp4 <- data.tmp3[data.tmp3$stateID == unique(region_master$stateID)[i],]
    clim.tmp[v] <- mean(as.vector(unlist(data.tmp4[,2:ncol(data.tmp4)])), na.rm = TRUE)
  }
  
  clim[[t]] <- clim.tmp
  
}

summary(lm(ar ~ unlist(ps) + unlist(lapply(clim, "[[", 1) ) + unlist(lapply(clim, "[[", 4))))






