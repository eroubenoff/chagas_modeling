library(lfe)
library(dplyr)
library(lme4)
library(ggplot2)
library(tidyr)
library(merTools)
library(RCurl)
library(pracma)


load("~/Google Drive/Dengue/Dengue-Malaria/Brazil/cleaned_data/yearly_all_munic.Rdata")

### time series of cases
ts <- yearly.dat.all %>% group_by(region, year) %>% summarize(urb = mean(urb, na.rm = T), 
                                                              gdp = mean(gdp, na.rm = T),
                                                              den = sum(den), ps = sum(ps, na.rm = T)) %>%
  mutate(den_inc = den*100000/ps) %>%
  arrange(region, year)

head(ts)


ggplot(ts, aes(x = year, y = den_inc, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  labs(x = "Year", y = "Yearly dengue incidence per 100,000",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())

ggplot(ts, aes(x = year, y = ps, group = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South"))
)) +
  geom_line(aes(color = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")), 
                linetype = factor(region, labels = c("North", "Northeast", "Central West", "Southeast", "South")))) +
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  labs(x = "Year", y = "",
       title = "", subtitle = "") + 
  theme(legend.position="top") +
  theme(legend.title=element_blank())



ts <- yearly.dat.all %>% group_by(stateID, year, region) %>% summarize(urb = mean(urb, na.rm = T), 
                                                              gdp = mean(gdp, na.rm = T),
                                                              den = sum(den), ps = sum(ps, na.rm = T)) %>%
  mutate(den_inc = den*100000/ps) %>%
  arrange(stateID, year)

head(ts)



ggplot(ts, aes(x = year, y = den_inc, group = factor(stateID))) +
  geom_line(aes(color = factor(stateID)), show.legend = FALSE)+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Spectral"))(colourCount)) +
  theme_minimal() +
  labs(x = "Year", y = "",
       title = "", subtitle = "") + 
  theme(legend.title=element_blank(), legend.position = "bottom") +
  facet_wrap(~region)


# simple model

model.mat <- yearly.dat.all
model.mat[,6:13] <- scale(yearly.dat.all[,6:13]) %>% as.data.frame()
model.mat <- model.mat %>% mutate(den_inc = den/ps, mal_inc = mal/ps)

formula <- den ~ offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI + year + factor(stateID)
formula.m <- mal ~ offset(log(exp.count.mal)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI + year + factor(stateID)


summary(glm(formula = formula.m, family = "poisson", data = model.mat))

# random effects model
form <- den ~  offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI  + (1 |code)
form.m <- mal ~  offset(log(exp.count.mal)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI  + (1 |code)

glm.fit <- glmer(form.m, data = model.mat, family = poisson(link = "log"))
summary(glm.fit)


nrow(model.mat)
length(getME(glm.fit,"theta"))
length(fixef(glm.fit))

tt <- getME(glm.fit,"theta")
ll <- getME(glm.fit,"lower")
min(tt[ll==0])

derivs1 <- glm.fit@optinfo$derivs
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))
max(pmin(abs(sc_grad1),abs(derivs1$gradient)))

dd <- update(glm.fit,devFunOnly=TRUE)
pars <- unlist(getME(glm.fit,c("theta","fixef")))
grad2 <- grad(dd,pars)
hess2 <- hessian(dd,pars)
sc_grad2 <- solve(hess2,grad2)
max(pmin(abs(sc_grad2),abs(grad2)))

ss <- getME(glm.fit,c("theta","fixef"))
m3 <- update(glm.fit,start=ss,control=glmerControl(optimizer="bobyqa",
                                                 optCtrl=list(maxfun=2e5)))

source(system.file("utils", "allFit.R", package="lme4"))

aa <- allFit(glm.fit)
is.OK <- sapply(aa,is,"merMod")  
aa.OK <- aa[is.OK]
lapply(aa.OK,function(x) x@optinfo$conv$lme4$messages)

(lliks <- sort(sapply(aa.OK,logLik)))




preds <- predictInterval(lme1, newdata = newDat, n.sims = 999)


form <- den ~ offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI | year + stateID | 0 | stateID
fit <- felm(form, data = model.mat)
summary(fit)

form <- mal ~ offset(log(exp.count.mal)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI| year:stateID | 0 | stateID
fit.mal <- felm(form, data = model.mat)
summary(fit.mal)

#form <- den ~ offset(log(exp.count)) + (1 + year|region)
#glm.fit <- glmer(form, data = model.mat, family = poisson(link = "log"))
#summary(glm.fit)


library(pglm)
fixed <- pglm(den ~ offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI, data=model.mat, 
              index=c("code", "year"), model="within",
              family = poisson(link = "log"))
summary(fixed)
              
random <-  pglm(den ~ offset(log(exp.count)) + tmp.av  + pre.av  + gdp + dens + urb + NDVI, data=model.mat, 
                index=c("code", "year"), model="random",
                family = poisson)
summary(random)



# library(plm)

# fixed <- plm(den_inc ~  tmp.av  + pre.av  + gdp + dens + urb + NDVI, 
#              data=model.mat, index=c("code", "year"), model="within")
# summary(fixed)
# 
# random <- plm(den_inc ~  tmp.av  + pre.av  + gdp + dens + urb + NDVI, 
#              data=model.mat, index=c("code", "year"), model="random")
# summary(random)
# 
# phtest(fixed, random)
# 
# 
# fixed.time <- plm(den_inc ~  tmp.av  + pre.av  + gdp + dens + urb + NDVI + factor(year), 
#                   data=model.mat, index=c("code", "year"), model="within")
# summary(fixed.time)
# pFtest(fixed.time, fixed)
# 
# fixed <- plm(mal_inc ~  tmp.av  + pre.av  + gdp + dens + urb + NDVI, 
#              data=model.mat, index=c("code", "year"), model="within")
# summary(fixed)
# 
# random <- plm(mal_inc ~  tmp.av  + pre.av  + gdp + dens + urb + NDVI, 
#               data=model.mat, index=c("code", "year"), model="random")
# summary(random)
# 
# phtest(fixed, random)
