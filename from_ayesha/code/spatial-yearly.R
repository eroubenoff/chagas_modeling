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

### collapse to yearly level
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic.Rdata")

r = tmp7 %>% group_by(year) %>% summarise(r = sum(den, na.rm = TRUE)/sum(ps, na.rm = TRUE))
r.mal = tmp7 %>% group_by(year) %>% summarise(r.mal = sum(mal, na.rm = TRUE)/sum(ps, na.rm = TRUE))

yearly.dat.all <- tmp7 %>% group_by(code, year) %>% summarize(den = sum(den,na.rm = T),
                                                              mal = sum(mal, na.rm = T),
                                                          ps = mean(ps, na.rm = T),
                                                          tmp.av = mean(tmp, na.rm = T),
                                                          tmp.max = max(tmp, na.rm = T),
                                                          pre.av = mean(pre, na.rm = T),
                                                          pre.max = max(pre, na.rm = T),
                                                          NDVI = mean(NDVI, na.rm = T),
                                                          gdp = mean(gdp, na.rm = T),
                                                          dens = mean(dens, na.rm = T),
                                                          urb = mean(urban.perc, na.rm = T),
                                                          stateID = unique(stateID),
                                                          region = unique(region))
                                                          
                                                          
                                                                  

yearly.dat.all <- left_join(yearly.dat.all, r, by = "year")
yearly.dat.all <- left_join(yearly.dat.all, r.mal, by = "year")

yearly.dat.all$exp.count <- yearly.dat.all$ps * yearly.dat.all$r
yearly.dat.all$exp.count.mal <- yearly.dat.all$ps * yearly.dat.all$r.mal

yearly.dat.all$mal[yearly.dat.all$year == 2007] <- NA


save(yearly.dat.all, file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/yearly_all_munic.Rdata")

#load shapefile
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/yearly_all_munic.Rdata")
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/")
municBRorig   <- readOGR(".", "municipios_2010")
municBRorig$code <- as.numeric(substr(municBRorig@data$codigo_ibg, 1,6))



# by region, specify time and region
unique(yearly.dat.all$region)

#merge south and southwest
yearly.dat.all$region[yearly.dat.all$region == 5] <- 4

for (r in 1:4){
  yearly.dat <- yearly.dat.all[yearly.dat.all$region == r, ]
  yearly.dat$year <- as.numeric(as.character(yearly.dat$year))
  
  yearly.dat <- yearly.dat[yearly.dat$year <= 2010, ] #no urbanization data after 2010
  
  municBR <- municBRorig
  municBR <- municBR[(municBR$code %in% yearly.dat$code),]
  municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$gdp))]),]
  municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$urb))]),]
  municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$NDVI))]),]
  municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$pre.av))]),]
  municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$tmp.av))]),]
  
  W.nb <- poly2nb(municBR, row.names = municBR$code)
  # drop places with no links
  sub_municBR <- subset(municBR, subset=card(W.nb) > 0)
  sub_W.nb <- subset(W.nb, subset=card(W.nb) > 0)
  W.list <- nb2listw(sub_W.nb, style = "B", zero.policy = TRUE)
  W <- nb2mat(sub_W.nb, style = "B", zero.policy = TRUE)
  
  # data subset
  sp.dat.sub <- subset(yearly.dat, yearly.dat$code %in% sub_municBR$code)
  sp.dat.sub <- sp.dat.sub[order(sp.dat.sub$year, sp.dat.sub$code),]
  head(sp.dat.sub)
  
  
  #standardize all covariates to zero mean and unit variance
  sp.dat.sub[,6:13] <- scale(sp.dat.sub[,6:13])
  
  
  #formula <- den ~ offset(log(exp.count)) + ns(tmp.av, df = 3)  + pre.av  + gdp + dens + urb + ns(NDVI, df = 3)
  formula <- den ~ offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI
  
  
  #summary(glm(formula = formula, family = "poisson", data = sp.dat.sub))
  
  model2 <- ST.CARar(formula = formula, family = "poisson",
                     data = sp.dat.sub, W = W, burnin = 20000, n.sample = 220000,
                     thin = 10)
  #plot(exp(model2$samples$beta[ ,1:3]))
  parameter.summary <- summarise.samples((model2$samples$beta[ , -1]),
                                         quantiles = c(0.5, 0.025, 0.975))
  df = data.frame(round(parameter.summary$quantiles, 5))
  
  save(df, file = paste0("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modelfit_region_", r, ".Rdata"))
  save(df, file = paste0("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modelfit_region_4+5", ".Rdata"))
  
  rm(model2)
  
}

#print(model2)
#summary(model2$samples)
#plot(exp(model2$samples$beta[ , -1]))
#parameter.summary <- summarise.samples((model2$samples$beta[ , -1]),quantiles = c(0.5, 0.025, 0.975))
# df = data.frame(round(parameter.summary$quantiles, 5))

region = c("North", "Northeast", "Central West", "South and Southeast")
df.all <- NA
for (r in 1:4){
  if (r == 4) {
    load(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modelfit_region_4+5.Rdata")
  } else {
    load(file = paste0("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/modelfit_region_", r, ".Rdata"))
    }
  
  names(df) <- c("est", "lower", "upper")
  df$est <- exp(df$est); df$lower <- exp(df$lower); df$upper <- exp(df$upper)
  #df$est <- (df$est); df$lower <- (df$lower); df$upper <- (df$upper)
  
  df$label <- c("Temp", "Prec", "GDP", "Density", "% Urban", "NDVI")
  df$region = region[r]
  df.all <- rbind(df.all, df)
}

df.all <- df.all[-1, ]
df.all$region <- factor(df.all$region, levels = c("North", "Northeast", "Central West", "South and Southeast"))

df.sub <- df.all[df.all$region != "South" & df.all$region != "Southeast", ]

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/model_results_by_region.pdf", width=8,height=6)

ggplot(data=df.all, aes(x=label, y=est, color = region)) +
  geom_point(size = 1, position=position_dodge(width=0.5)) + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5, position=position_dodge(width=0.5))+
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Relative risk") +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())

dev.off()

########################################################
# malaria

load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/yearly_all_munic.Rdata")


yearly.dat <- yearly.dat.all[yearly.dat.all$region %in% r, ]
yearly.dat$year <- as.numeric(as.character(yearly.dat$year))

yearly.dat <- yearly.dat[yearly.dat$year <= 2010, ] #no urbanization data after 2010

municBR <- municBRorig
municBR <- municBR[(municBR$code %in% yearly.dat$code),]
municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$gdp))]),]
municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$urb))]),]
municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$NDVI))]),]
municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$pre.av))]),]
municBR <- municBR[!(municBR$code %in% yearly.dat$code[which(is.na(yearly.dat$tmp.av))]),]


W.nb <- poly2nb(municBR, row.names = municBR$code)
# drop places with no links
sub_municBR <- subset(municBR, subset=card(W.nb) > 0)
sub_W.nb <- subset(W.nb, subset=card(W.nb) > 0)
W.list <- nb2listw(sub_W.nb, style = "B", zero.policy = TRUE)
W <- nb2mat(sub_W.nb, style = "B", zero.policy = TRUE)

# data subset
sp.dat.sub <- subset(yearly.dat, yearly.dat$code %in% sub_municBR$code)
sp.dat.sub <- sp.dat.sub[order(sp.dat.sub$year, sp.dat.sub$code),]

#standardize all covariates to zero mean and unit variance
sp.dat.sub[,6:13] <- scale(sp.dat.sub[,6:13])


formula <- mal ~ offset(log(exp.count.mal)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI



summary(glm(formula = formula, family = "poisson",data = sp.dat.sub))




model2 <- ST.CARar(formula = formula, family = "poisson",
                   data = sp.dat.sub, W = W, burnin = 20000, n.sample = 220000,
                   thin = 10)
print(model2)
summary(model2$samples)
plot(exp(model2$samples$beta[ , -1]))

parameter.summary <- summarise.samples((model2$samples$beta[ , -1]),
                                       quantiles = c(0.5, 0.025, 0.975))
