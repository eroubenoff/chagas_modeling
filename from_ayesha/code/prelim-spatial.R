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

rm(list = ls())

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic.Rdata")
sp.dat <- tmp7[,c("den", "exp.count", "month", "year", "tmp", "pre", "NDVI", "gdp","dens", "lag1_tmp", "lag1_pre", "stateID", "region", "code") ]





sp.dat$SMR <- sp.dat$den/sp.dat$exp.count
sp.dat$logSMR <- log(sp.dat$SMR)
sp.dat <- sp.dat[order(sp.dat$year, sp.dat$month, sp.dat$code),]

#par(pty="s", cex.axis=1.5, cex.lab=1.5)
#pairs(sp.dat[ , c("logSMR", "tmp", "pre", "NDVI")], pch = 19, cex = 0.5,
#         lower.panel=NULL, panel=panel.smooth,
#         labels = c("ln(SMR)", "temp", "prec", "NDVI"))

#load shapefile
municBR   <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/municipios_2010.shp")
statesBR  <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/estados_2010/estados_2010.shp")
regionsBR <- readShapePoly(fn="~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/regioes_2010/regioes_2010.shp")
municBR@data$code <- as.numeric(substr(municBR@data$codigo_ibg, 1,6))
municBR@data$order <- 1:nrow(municBR@data)

SMR.av <- sp.dat %>% group_by(code) %>% summarize(SMR.mean = mean(SMR))
municBR@data <- left_join(municBR@data, SMR.av, by = "code")
 

# l1 = list("SpatialPolygonsRescale", layout.north.arrow(),
#              offset = c(220000,647000), scale = 4000)
# l2 = list("SpatialPolygonsRescale", layout.scale.bar(), offset =
#           c(225000, 647000), scale = 10000, fill = c("transparent","black"))
# l3 = list("sp.text", c(225000,649000), "0")
# l4 = list("sp.text", c(230000,649000), "5000 m")
# breakpoints <- seq(min(SMR.av$SMR.mean, na.rm = T)-0.1, max(SMR.av$SMR.mean, na.rm = T)+0.1,
#                   length.out = 11)
# spplot(municBR, "SMR.mean", sp.layout = list(l1, l2, l3, l4),
#        xlab = "Easting", ylab = "Northing",
#        scales = list(draw = TRUE),  at = breakpoints,
#        col.regions = terrain.colors(n = length(breakpoints)-1),
#        par.settings=list(fontsize=list(text=20)))

municBR <- municBR[municBR$code != 250860, ] # missing tmp and pre for this code
W.nb <- poly2nb(municBR, row.names = municBR$code)
#W.list <- nb2listw(W.nb, style = "B", zero.policy = TRUE)
#W <- nb2mat(W.nb, style = "B", zero.policy = TRUE)

# drop places with no links
sub_municBR <- subset(municBR, subset=card(W.nb) > 0)
sub_W.nb <- subset(W.nb, subset=card(W.nb) > 0)
W.list <- nb2listw(sub_W.nb, style = "B", zero.policy = TRUE)
W <- nb2mat(sub_W.nb, style = "B", zero.policy = TRUE)

# data subset
sp.dat.sub <- subset(sp.dat, sp.dat$code %in% sub_municBR$code)
sp.dat.sub$year <- as.numeric(as.character(sp.dat.sub$year))
sp.dat.sub <- sp.dat.sub[sp.dat.sub$year > 2005, ]

which(is.na(sp.dat.sub$gdp))
sp.dat.sub$year[which(is.na(sp.dat.sub$gdp))]
sp.dat.sub[which(is.na(sp.dat.sub$exp.count)), ][1:5,]


formula <- den ~ offset(log(exp.count)) + tmp + pre + gdp + dens

model1 <- glm(formula = formula, family = "quasipoisson",
                 data = sp.dat.sub)
resid.glm <- residuals(model1)
summary(model1)$coefficients

summary(model1)$dispersion

summary(model1)

model2 <- ST.CARar(formula = formula, family = "poisson",
                      data = sp.dat.sub, W = W, burnin = 1000, n.sample = 10000,
                      thin = 10)
print(model2)
summary(model2$samples)
plot(exp(model2$samples$beta[ , -1]))
