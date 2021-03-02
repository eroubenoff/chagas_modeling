library("sp")
library(dplyr)
library("spdep")
library(CARBayesST)
library(tidyr)
library(mgcv)
library(splines)
library(ggplot2)
library(rgdal)
library(maptools)
require(descr)
require(RColorBrewer)
require(plotrix)
library(lfe)
library("CARBayes")

rm(list = ls())

### create state level file
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic.Rdata")
r = tmp7 %>% group_by(month) %>% summarise(r = sum(den, na.rm = TRUE)/sum(ps, na.rm = TRUE))


state.dat <- tmp7 %>% group_by(stateID, year, month) %>% summarize(den = sum(den,na.rm = T),
                                                                   ps = sum(ps, na.rm = T),
                                                                   tmp = mean(tmp, na.rm = T),
                                                                   pre = mean(pre, na.rm = T),
                                                                   NDVI = mean(NDVI, na.rm = T),
                                                                   gdp = mean(gdp, na.rm = T),
                                                                   dens = mean(dens, na.rm = T),
                                                                   region = unique(region))

state.dat <- state.dat %>% group_by(stateID) %>% mutate(lag1_pre = lag(pre), lag1_tmp = lag(tmp))
state.dat <- left_join(state.dat, r, by = "month")
state.dat$exp.count <- state.dat$ps * state.dat$r

sp.dat <- state.dat
sp.dat$SMR <- sp.dat$den/sp.dat$exp.count
sp.dat$logSMR <- log(sp.dat$SMR)
sp.dat <- sp.dat[order(sp.dat$year, sp.dat$month, sp.dat$stateID),]

#load shapefile
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
statesBR@data$stateID <- as.numeric(as.character(statesBR@data$codigo_ibg))
statesBR@data$order <- 1:nrow(statesBR@data)

SMR.av <- sp.dat %>% group_by(stateID) %>% summarize(SMR.mean = mean(SMR))
SMR.av$stateID <- as.numeric(as.character(SMR.av$stateID))
statesBR@data <- left_join(statesBR@data, SMR.av, by = "stateID")

breakpoints <- seq(min(SMR.av$SMR.mean, na.rm = T)-0.1, max(SMR.av$SMR.mean, na.rm = T)+0.1,
                   length.out = 11)

spplot(statesBR, "SMR.mean",
        scales = list(draw = TRUE),  at = breakpoints,
        col.regions = terrain.colors(n = length(breakpoints)-1),
        par.settings=list(fontsize=list(text=10)))

W.nb <- poly2nb(statesBR, row.names = statesBR$stateID)
W.list <- nb2listw(W.nb, style = "B", zero.policy = TRUE)
W <- nb2mat(W.nb, style = "B", zero.policy = TRUE)

# data subset
sp.dat$year <- as.numeric(as.character(sp.dat$year))
sp.dat.sub <- sp.dat[sp.dat$year > 2007, ]



formula <- den ~ offset(log(exp.count)) + tmp + pre + gdp + dens

model1 <- glm(formula = formula, family = "quasipoisson",
              data = sp.dat.sub)
resid.glm <- residuals(model1)
summary(model1)$coefficients

summary(model1)$dispersion

summary(model1)

model2 <- ST.CARar(formula = formula, family = "poisson",
                   data = sp.dat.sub, W = W, burnin = 20000, n.sample = 220000,
                   thin = 10)
print(model2)
summary(model2$samples)
plot(exp(model2$samples$beta[ , -1]))
parameter.summary <- summarise.samples(exp(model2$samples$beta[ , -1]),
                                       quantiles = c(0.5, 0.025, 0.975))
round(parameter.summary$quantiles, 3)

