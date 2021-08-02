library(Rcpp)
library(gridExtra)
library(ggplot2)
library(lme4)
library(dplyr)
library(survival)
library(survminer)
library(Hmisc)
library(PerformanceAnalytics)
library(gnlm)
library(MuMIn)
library(LMERConvenienceFunctions)
library(lsmeans)
library(effects)
library(olsrr)

kites <- read.csv("kite_data.csv", header = T, fileEncoding = "UTF-8-BOM")

head(kites)
glimpse(kites)

#standardise numeric explanatory varibles
X._urban_HR2 <- scale(kites$X._urban_HR, center = TRUE, scale = TRUE)
X._urban_core2 <- scale(kites$X._urban_core, center = TRUE, scale = TRUE)
X._closed_HR2 <- scale(kites$X._closed_HR, center = TRUE, scale = TRUE)
X._closed_core2 <- scale(kites$X._closed_core, center = TRUE, scale = TRUE)
X._open_HR2 <- scale(kites$X._open_HR, center = TRUE, scale = TRUE)
X._open_core2 <- scale(kites$X._open_core, center = TRUE, scale = TRUE)
elev_range_HR <- scale(kites$elev_range_HR_km, center = TRUE, scale = TRUE)
elev_range_core <- scale(kites$elev_range_core_km, center = TRUE, scale = TRUE)
prop_low_HR <- scale(kites$prop_low_HR, center = TRUE, scale = TRUE)
prop_low_core <-scale(kites$prop_low_core, center = TRUE, scale = TRUE)
winter_duration_days2 <- scale(kites$winter_duration_days, center = TRUE, scale = TRUE)
gps_fixes2 <- scale(kites$gps_fixes, center = TRUE, scale = TRUE)

#Correlation analyses
#Distance travelled (variables are the same for home range (95%))
dist_contin_kites <- kites[, c(7,8,9,10,11,12,15,16,17,18,21,22)] #isolate the continuous explanatory variables
head(dist_contin_kites) #double check they are all included

res1 <- rcorr(as.matrix(dist_contin_kites)) #correlation matrix command
res1 #show results

res1$r #extract the correlation coefficients
res1$P #extract the p-values

mytab <- res1$r
write.table(mytab, file = "distance_correlation_coefficients.csv")

chart.Correlation(dist_contin_kites, histrogram=TRUE, pch=22, cex.cor.scale=2)

#Core areas (50% KDE)
core_contin_kites <- kites[, c(8,10,12,14,17,18)] #isolate the continuous explanatory variables
head(core_contin_kites) #double check they are all included

res2 <- rcorr(as.matrix(core_contin_kites)) #correlation matrix command
res2 #show results

res2$r #extract the correlation coefficients
res2$P #extract the p-values

mytab <- res2$r
write.table(mytab2, file = "core_areas_correlation_coefficients.csv")

chart.Correlation(core_contin_kites, histrogram=TRUE, pch=19)

#Remove outliers (see manuscript main text 'Methods' for details)
kites_noouts <- read.csv("kite_data_outliers_removed.csv", header = T, fileEncoding = "UTF-8-BOM")

#rescale continuous explanatory variables
X._urban_HR2 <- scale(kites$X._urban_HR, center = TRUE, scale = TRUE)
X._open_HR2 <- scale(kites$X._open_HR, center = TRUE, scale = TRUE)
prop_low_HR2 <- scale(kites$prop_low_HR, center = TRUE, scale = TRUE)
winter_duration_days2 <- scale(kites$winter_duration_days, center = TRUE, scale = TRUE)

X._urban_HR3 <- scale(kites_noouts$X._urban_HR, center = TRUE, scale = TRUE)
X._open_HR3 <- scale(kites_noouts$X._open_HR, center = TRUE, scale = TRUE)
prop_low_HR3 <- scale(kites_noouts$prop_low_HR, center = TRUE, scale = TRUE)
winter_duration_days3 <- scale(kites_noouts$winter_duration_days, center = TRUE, scale = TRUE)

X._urban_core4 <- scale(kites_noouts$X._urban_core, center = TRUE, scale = TRUE)
X._open_core4 <- scale(kites_noouts$X._open_core, center = TRUE, scale = TRUE)
prop_low_core4 <- scale(kites_noouts$prop_low_core, center = TRUE, scale = TRUE)
winter_duration_days4 <- scale(kites_noouts$winter_duration_days, center = TRUE, scale = TRUE)

#AIC Model selection for standardized distance travelled models

mod0 <- glmer(std_dist_km ~ + (1|bird_id), data = kites, family = Gamma(link = log), nAGQ=0)

summary(mod0)

mod1 <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region + arrival_date + winter_duration_days2
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod1)

mod2 <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region + arrival_date
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod2)

mod3 <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod3)

mod4 <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2 + prop_low_HR2
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod4)

mod5 <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod5)

mod6 <- glmer(std_dist_km ~ sex + age + X._urban_HR2
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod6)

mod7 <- glmer(std_dist_km ~ sex + age
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod7)

mod8 <- glmer(std_dist_km ~ sex
              + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod8)

model.sel(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) #used this for model selection
dist_output <- model.sel(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8)
anova(mod0, mod1, mod2, mod3, mod4, mod5, mod6, mod7, mod8) #did not use this
write.table(dist_output, file = "std_model_selection.csv")

#Distance model coefficients and confidence intervals
MA.ests1 <- model.avg(dist_output)
MA.ests1
coefTable(MA.ests1)
confint(MA.ests1)
std_dist_coef <- cbind(coefTable(MA.ests1), confint(MA.ests1))
write.table(std_dist_coef, file = "std_dist_coef.csv")

#AIC Model selection for home range (95% KDE)

mod00 <- glmer(X95._HR_km2 ~ 
                 + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod00)

mod9 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region + arrival_date + winter_duration_days3
              + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod9)

mod10 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region + arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod10)

mod11 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod11)

mod12 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3 + prop_low_HR3
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod12)

mod13 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod13)

mod14 <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod14)

mod15 <- glmer(X95._HR_km2 ~ sex + age
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod15)

mod16 <- glmer(X95._HR_km2 ~ sex
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod16)

model.sel(mod00, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16)
anova(mod00, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16)
home_output <- model.sel(mod00, mod9, mod10, mod11, mod12, mod13, mod14, mod15, mod16)

MA.ests2 <- model.avg(home_output)
MA.ests2
coefTable(MA.ests2)
confint(MA.ests2)
cbind(coefTable(MA.ests2), confint(MA.ests2))

#AIC Model selection for core areas (50% KDE)

mod000 <- glmer(X50._core_km2 ~ 
                  + (1|bird_id), data = kites_noouts, family = Gamma(link=log))
summary(mod000)

mod17 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4 + prop_low_core4 + region + arrival_date + winter_duration_days4
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod17)

mod18 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4 + prop_low_core4 + region + arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod18)

mod19 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4 + prop_low_core4 + region
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod19)

mod20 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4 + prop_low_core4
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod20)

mod21 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod21)

mod22 <- glmer(X50._core_km2 ~ sex + age + X._urban_core4
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod22)

mod23 <- glmer(X50._core_km2 ~ sex + age
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod23)

mod24 <- glmer(X50._core_km2 ~ sex
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod24)

model.sel(mod000, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24)
anova(mod000, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24)
core_output <- model.sel(mod000, mod17, mod18, mod19, mod20, mod21, mod22, mod23, mod24)

MA.ests3 <- model.avg(core_output)
MA.ests3
coefTable(MA.ests3)
confint(MA.ests3)
cbind(coefTable(MA.ests3), confint(MA.ests3))

#Final generalized linear-mixed models

mdistance <- glmer(std_dist_km ~ sex + age + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region + arrival_date + winter_duration_days2
                   + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)

summary(mdistance)
anova(mdistance)
r.squaredGLMM(mdistance)
ols_test_f(mdistance)



mhome <- glmer(X95._HR_km2 ~ sex + age + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region + arrival_date + winter_duration_days3
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)

summary(mhome)
anova(mhome)
r.squaredGLMM(mhome)
ols_test_f(mhome)



mcore <- glmer(X50._core_km2 ~ sex + age + X._urban_core4 + X._open_core4 + prop_low_core4 + region
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mcore)
anova(mcore)
r.squaredGLMM(mcore)
ols_test_f(mcore)
lsmeans(mcore, pairwise ~ region)

#Effects plots
plot(allEffects(mdistance), lty = 2, lwd = 2, multiline = T, main = NA, cex =5, colors = palette(), z.var = T)
plot(allEffects(mhome), lty = 2, lwd = 2, multiline = T, main = NA, cex =5, colors = palette(), z.var = T)
plot(allEffects(mcore), lty = 2, lwd = 2, multiline = T, main = NA, cex =5, colors = palette(), z.var = T)

#Interaction models
#AIC Model selection for standardized distance travelled models

mod0000 <- glmer(std_dist_km ~ + (1|bird_id), data = kites, family = Gamma(link = log), nAGQ=0)

summary(mod0000)

mod34 <- glmer(std_dist_km ~ sex + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region + winter_duration_days2 + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod34)
anova(mod34)
lsmeans(mod34, pairwise ~ age*arrival_date)

mod35 <- glmer(std_dist_km ~ sex + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + region + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod35)

mod36 <- glmer(std_dist_km ~ sex + X._urban_HR2 + X._open_HR2 + prop_low_HR2 + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod36)

mod37 <- glmer(std_dist_km ~ sex + X._urban_HR2 + X._open_HR2 + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod37)

mod38 <- glmer(std_dist_km ~ sex + X._urban_HR2 + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod38)

mod39 <- glmer(std_dist_km ~ sex + age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod39)

mod40 <- glmer(std_dist_km ~ sex
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
summary(mod40)

model.sel(mod0000, mod34, mod35, mod36, mod37, mod38, mod39, mod40)
anova(mod0000, mod34, mod35, mod36, mod37, mod38, mod39, mod40)
dist_interaction <- model.sel(mod0000, mod34, mod35, mod36, mod37, mod38, mod39, mod40)
write.table(dist_interaction, file = "distance_interaction_modsel.csv")


MA.ests5 <- model.avg(dist_interaction)
MA.ests5
coefTable(MA.ests5)
confint(MA.ests5)
dist_int_coefficients <- cbind(coefTable(MA.ests5), confint(MA.ests5))
write.table(dist_int_coefficients, file = "dist_int_coef.csv")

mod41 <- glmer(std_dist_km ~ age*arrival_date
               + (1|bird_id), data = kites, family = Gamma(link=log), nAGQ=0)
plot(predictorEffects(mod41, ~ arrival_date, partial.residuals=T))

#AIC Model selection for home range (95% KDE) interaction models

mod00000 <- glmer(X95._HR_km2 ~ 
                    + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod00000)

mod42 <- glmer(X95._HR_km2 ~ sex + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region + winter_duration_days3 + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod42)

mod43 <- glmer(X95._HR_km2 ~ sex + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + region + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod43)

mod44 <- glmer(X95._HR_km2 ~ sex + X._urban_HR3 + X._open_HR3 + prop_low_HR3 + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod44)

mod45 <- glmer(X95._HR_km2 ~ sex + X._urban_HR3 + X._open_HR3 + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod45)

mod46 <- glmer(X95._HR_km2 ~ sex + X._urban_HR3 + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod46)

mod47 <- glmer(X95._HR_km2 ~ sex + age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod47)

mod48 <- glmer(X95._HR_km2 ~ sex
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
summary(mod48)

model.sel(mod00000, mod42, mod43, mod44, mod45, mod46, mod47, mod48)
anova(mod00000, mod42, mod43, mod44, mod45, mod46, mod47, mod48)
home_range_int <- model.sel(mod00000, mod42, mod43, mod44, mod45, mod46, mod47, mod48)
write.table(home_range_int, file = "home_range_int_modelsel.csv")


MA.ests6 <- model.avg(home_range_int)
MA.ests6
coefTable(MA.ests6)
confint(MA.ests6)
home_range_int_coefficients <- cbind(coefTable(MA.ests6), confint(MA.ests6))
write.table(home_range_int_coefficients, file = "home_range_int_coef.csv")

#Final home range (95% KDE) interaction model
mod49 <- glmer(X95._HR_km2 ~ age*arrival_date
               + (1|bird_id), data = kites_noouts, family = Gamma(link=log), nAGQ=0)
plot(predictorEffects(mod49, ~ arrival_date, partial.residuals=T))

#Survival analyses
#Fisher's exact test on kite survival
surv_knowns <- read.csv("kite_survival_data.csv")

sex_tab <- table(surv_knowns$sex, surv_knowns$death_type)
barplot(sex_tab, beside = T, legend = T)
sex_fishers <- fisher.test(sex_tab, conf.int = T, conf.level = 0.95)
sex_fishers

region_tab <- table(surv_knowns$winter_region, surv_knowns$death_type)
barplot(region_tab, beside = T, legend = T)
region_fishers <- fisher.test(region_tab, conf.int = T, conf.level = 0.95)
region_fishers

