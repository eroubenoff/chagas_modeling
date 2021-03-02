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
library(tmaptools)
library(tmap)

rm(list = ls())

### collapse to yearly level
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/tmp_all_munic.Rdata")
load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/yearly_all_munic.Rdata")
####

### time series of cases
ts <- tmp7 %>% group_by(region, year, month) %>% summarize(den = sum(den), ps = sum(ps, na.rm = T)) %>%
  mutate(den_inc = den*100000/ps)

ts <- ts[order(ts$region, ts$year, ts$month), ]
head(ts)

ts$time <- rep(seq(2001, 2013, by = 1/12)[-length(seq(2001, 2013, by = 1/12))], 5)

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/timeseries_by_region.pdf", width=8,height=6)
ggplot(ts, aes(x = time, y = den_inc, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
               )) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_classic() +
  labs(x = "Month", y = "Monthly dengue incidence per 100,000",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())
dev.off()

# NDVI
yearly.dat.all %>% group_by(region) %>% summarize(n = length(unique(stateID)))

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_NDVI_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=NDVI, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "NDVI", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
dev.off()

states <- tmp7 %>% group_by(stateID, stateName) %>% summarize(n())
ggplot(yearly.dat.all[yearly.dat.all$region == "3" ,], aes(x=NDVI, y=den*100000/ps, col = stateID)) +
  geom_point(alpha = 0.5)  


# urbanization
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_urban_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=urb, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "% urban", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()

# average temp for the year
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_temp_avg_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=tmp.av, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "Average temp (Celsius)", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()

# average prec for the year
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_prec_avg_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=pre.av, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "Average precipitation", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()

# GDP
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_gdp_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=gdp, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "GDP", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()

# density
png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/den_dens_reg.png", width=5,height=8,units="in",res=1200)
ggplot(yearly.dat.all, aes(x=dens, y=den*100000/ps)) +
  geom_point(alpha = 1/5) + 
  theme(axis.text=element_text(size=7),axis.title=element_text(size=7), strip.text.x = element_text(size=7, margin = margin(.08, 0, .08, 0, "cm"))) +
  facet_wrap(~ region, ncol = 1)+ labs(x = "Density", y = "dengue cases per 100,000") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

dev.off()


#plot yearly incidence
yearly.dat.all$den_inc <- yearly.dat.all$den*100000/yearly.dat.all$ps
setwd("~/Google Drive/Dengue/Dengue-Malaria/Brazil/mapdata/municipios_2010/")
municBRorig   <- readOGR(".", "municipios_2010")
municBRorig$code <- as.numeric(substr(municBRorig@data$codigo_ibg, 1,6))


inc.av <- yearly.dat.all %>% group_by(code) %>% summarize(inc.mean = mean(den_inc, na.rm = TRUE))
inc.max <- yearly.dat.all %>% group_by(code) %>% summarize(inc.max = max(den_inc, na.rm = TRUE))
municBR <- municBRorig
municBR@data <- left_join(municBR@data, inc.av, by = "code")
municBR@data <- left_join(municBR@data, inc.max, by = "code")


municBR@data$inc.mean[municBR@data$inc.mean == 0] <- NA

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_dengue_inc.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("inc.mean", textNA = "No cases", colorNA = "grey", 
          breaks = c(1,250,500,1000,1500,2000,2500), palette = "YlOrRd", title= "Average yearly dengue incidence per 100,000") 
dev.off()

municBR@data$inc.max[municBR@data$inc.max == 0] <- NA

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_dengue_inc_max.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("inc.max", textNA = "No cases", colorNA = "grey", 
          breaks = c(1,500,1000,2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000),palette = "YlOrRd", title= "Max dengue incidence per 100,000 (2001-2010)") 
dev.off()


#plot yearly avg temp
year.av <- yearly.dat.all %>% group_by(code) %>% summarize(year.av = mean(tmp.av, na.rm = TRUE))
municBR <- municBRorig
municBR@data <- left_join(municBR@data, year.av, by = "code")

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_temp.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("year.av", 
           palette = "YlOrRd", title= "Average yearly temperature") 
dev.off()

#plot yearly avg prec
year.av <- yearly.dat.all %>% group_by(code) %>% summarize(year.av = mean(pre.av, na.rm = TRUE))
municBR <- municBRorig
municBR@data <- left_join(municBR@data, year.av, by = "code")

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_prec.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("year.av", 
          palette = "Blues", title= "Average yearly precipitation") 
dev.off()


#plot yearly % urban
year.av <- yearly.dat.all %>% group_by(code) %>% summarize(year.av = mean(urb, na.rm = TRUE))
municBR <- municBRorig
municBR@data <- left_join(municBR@data, year.av, by = "code")

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_urb.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("year.av", style = "quantile",
          palette = "Blues", title= "% Urban (average across 10 years)") 
dev.off()

#plot avg pop
year.av <- yearly.dat.all %>% group_by(code) %>% summarize(year.av = mean(dens, na.rm = TRUE))
municBR <- municBRorig
municBR@data <- left_join(municBR@data, year.av, by = "code")
summary(municBR@data$year.av)

png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/map_density.png", width=6,height=8,units="in",res=1200)
tm_shape(municBR) +
  tm_fill("year.av", style = "quantile",
          palette = "Blues", title= "Average density per km2 from 2001-2010 (in quantiles)") 
dev.off()


### seasonality by month
monthly_inc <- tmp7 %>% group_by(region, month) %>% summarise(mean = mean(den_inc, na.rm = TRUE), 
                                                              sd = sd(den_inc, na.rm = TRUE)) %>% 
  group_by(region) %>% mutate(mean.reg = mean(mean)) %>% ungroup() %>% mutate(mean.centered = mean - mean.reg) 

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/seasonality_dengue_region.pdf", width=8,height=7)
ggplot(monthly_inc, aes(x=as.numeric(month), y=mean, color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.1) +
  geom_line() +
  #scale_colour_hue(l = 60) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() +
  scale_x_continuous(breaks=c(1:12)) +
  labs(x = "Month", y = "Mean dengue incidence per 100,000", color = "Region",
       title = "", subtitle = "") +
  theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) 



dev.off()


png(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/seasonality_meanCentered_dengue_region.png", width=8,height=7,units="in",res=1200)
ggplot(monthly_inc, aes(x=as.numeric(month), y=mean.centered, color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.1) +
  geom_line() +
  #scale_colour_hue(l = 60) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() +
  scale_x_continuous(breaks=c(1:12)) +
  labs(x = "Month", y = "Mean dengue incidence per 100,000", color = "Region",
       title = "", subtitle = "") +
  theme(legend.position="top") +
  theme(legend.title=element_blank())

dev.off()



monthly_inc <- tmp7 %>% group_by(region, month) %>% summarise(mean = mean(tmp, na.rm = TRUE), 
                                                              sd = sd(tmp, na.rm = TRUE)) %>% 
  group_by(region) %>% mutate(mean.reg = mean(mean)) %>% ungroup() %>% mutate(mean.centered = mean - mean.reg) 

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/seasonality_temp_region.pdf", width=8,height=7)
ggplot(monthly_inc, aes(x=as.numeric(month), y=mean, color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.1) +
  geom_line() +
  geom_hline(yintercept = 26, linetype = 2) + 
  geom_hline(yintercept = 18, linetype = 2) +
  annotate("text",10, 25, vjust = -1, label = "Peak Ae. aegypti transmission") +
  #scale_colour_hue(l = 60) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() +
  scale_x_continuous(breaks=c(1:12)) +
  labs(x = "Month", y = "Mean temperature (Celsius)", color = "Region",
       title = "", subtitle = "") +
  theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) 

# most transmission by Ae. aegypti
# transmission occurs between 18-34, peak at 26
dev.off()



monthly_inc <- tmp7 %>% group_by(region, month) %>% summarise(mean = mean(pre, na.rm = TRUE), 
                                                              sd = sd(pre, na.rm = TRUE)) %>% 
  group_by(region) %>% mutate(mean.reg = mean(mean)) %>% ungroup() %>% mutate(mean.centered = mean - mean.reg) 

pdf(file = "~/Google Drive/Dengue/Dengue-Malaria/Brazil/results/plots/seasonality_prec_region.pdf", width=8,height=7)
ggplot(monthly_inc, aes(x=as.numeric(month), y=mean, color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  geom_point() +
  #geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width=.1) +
  geom_line() +
  #scale_colour_hue(l = 60) +
  scale_colour_brewer(palette = "Dark2") +
  theme_classic() +
  scale_x_continuous(breaks=c(1:12)) +
  labs(x = "Month", y = "Mean precipitation", color = "Region",
       title = "", subtitle = "") +
  theme(legend.position="top", legend.title=element_blank(), legend.text=element_text(size=20), text = element_text(size=20)) 

dev.off()


