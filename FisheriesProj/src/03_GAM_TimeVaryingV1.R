library(brms)
library(ggplot2)
library(dplyr)
library(MARSS)
library(ggplot2)
library(broom.mixed)
library(MCMCvis)
library(tidyverse)
library(dplyr)
library(mgcv)
library(MASS)
library(stringr)
library(gamm4)
library(tidyr)
library(ggthemes)


SalmonSurv <- read.csv("data/Salmon/Survival_combined.csv")
SalmonSurvCUI<- SalmonSurv%>% #grabbing a single stock, I chose at random
  group_by(stock_name)%>%
  mutate(Survival_scale=log(Marine.Survival))%>% #making sure data are normalized
  ungroup()%>%
  rename(year=brood_year)%>%
  as.data.frame()%>%
  filter(calendar_year>=1986 & calendar_year<=2022&stock_name!='NA')
unique(SalmonSurvCUI$stock_name)
max_tree<-13
adapt_d <- 0.99


colnames(SalmonSurvCUI)
column.names<- c("Beuti_tumi_scale","pdo_spring", "sst_spring")
conditions=data.frame(calendar_year=unique(SalmonSurvCUI$calendar_year))
#conditions=data.frame(calendar_year=2013)

conditional_eff<-NA
bayes_r2<-NA
splines<-NA
stocks<-unique(na.omit(data.frame(SalmonSurvCUI))$stock_name)
for(j in 1:length(column.names)){
  cov_select<-  na.omit(data.frame(SalmonSurvCUI))
names(cov_select)[names(cov_select) == column.names[j]] <- 'covariate'
for(i in 1:length(stocks)){#length(stocks)
temp_surv<- cov_select%>%filter(stock_name==stocks[i])
initial_fit <- brm(Survival_scale ~  s(calendar_year, k=8) + #s(year,m=2, k=8) +
                     s(covariate, k=8, by =calendar_year),#+arma(year, p = 1),
                   data = temp_surv,
                   control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
                   iter = 3000, cores=3,chains = 3)

#effect for each year
temp_cond_eff<-conditional_effects(initial_fit, conditions = conditions)$covariate%>%
  mutate(stock_name=stocks[i], covariate_name =column.names[j])
conditional_eff<-rbind(temp_cond_eff,conditional_eff)
#bayes R2 value (fit)
temp_bayes_R2<-data.frame(bayes_R2(initial_fit))%>%
  mutate(stock_name=stocks[i], covariate_name =column.names[j])
bayes_r2<-rbind(temp_bayes_R2,bayes_r2)
#Spline SDs (linear or not)
temp_splines<-data.frame(summary(initial_fit)$splines)%>%
  mutate(stock_name=stocks[i],  covariate_name =column.names[j], spline=c('year', 'covariate'))
splines<-rbind(temp_splines,splines)
print(i)
print(column.names[j])
}
}

conditional_effects(initial_fit)$year

Identifiers<-SalmonSurvCUI%>%
  dplyr::select(RunType, stock_code, stock_name, area, ecoregion, Code, Code2)%>%
  distinct()
unique(SalmonSurvCUI$stock_name)
conditional_eff_id<-left_join(conditional_eff,Identifiers)%>%
  left_join(Atlas)
splines_id<-left_join(na.omit(splines),Identifiers)%>%
  left_join(Atlas)%>%add_column(spline=rep(c('year', 'covariate'), length(splines$Estimate)/2))
bayes_r2_id<-left_join(na.omit(bayes_r2),Identifiers)%>%
  left_join(Atlas)

write.csv(splines_id, file = "data/Salmon/splines_ALtCov1.csv")
write.csv(conditional_eff_id, file = "data/Salmon/conditional_eff_ALtCov1.csv")
write.csv(bayes_r2_id, file = "data/Salmon/bayes_r2_ALtCov1.csv")


filter(conditional_eff, covariate <= 0.03 & covariate >= -0.03)
conditional_eff<-conditional_effects(initial_fit, conditions = conditions)[column.names[j]]
plot(marginal_effects(initial_fit))
unique(conditional_eff$stock_name)

ggplot(filter(conditional_eff_id, covariate <= 0.03 & covariate >= -0.03),
       aes(x=calendar_year, y=estimate__, group=stock_name)) +
  theme_bw() +
  facet_grid(region~covariate_name)+
    geom_line(col='grey')

#### Spline plots ####

ggplot(splines_id,aes(x=RunType, y=Estimate, fill=RunType)) +
  theme_bw() +
  facet_grid(spline~covariate_name)+
    geom_point()+
  #scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) 

ggplot(splines_id,aes(x=region, y=Estimate, fill=region)) +
  theme_bw() +
  facet_grid(spline~covariate_name)+
    geom_point()+
  #scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) 

#### R2 plots ####

ggplot(bayes_r2_id,aes(x=RunType, y=Estimate, fill=RunType)) +
  theme_bw() +
  facet_wrap(~covariate_name)+
    geom_point()+
  #scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) 

ggplot(bayes_r2_id,aes(x=region, y=Estimate, fill=region)) +
  theme_bw() +
  facet_wrap(~covariate_name)+
  geom_point()+
  #scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) 




colnames(SalmonSurvCUI)
unique(SalmonSurvCUI$stock_name)
unique(SalmonSurvCUI$area)


# evaluate fit
pred <- predict(initial_fit)
joined <- cbind(SalmonSurvCUI, pred)
ggplot(joined, aes(year, Estimate)) + 
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),alpha=0.2) + 
  geom_line() + 
  geom_point(aes(year,Survival_scale)) + 
  ggtitle("Overall fit")

# To do - add this to the loop once you get it to run and practice plots for all 
# covariates
SalmonSurvCUI$predictions <- NA
nT <- nrow(SalmonSurvCUI)
for (i in (nT-10+1):nT) {
  # Update the model to use data up to the current time point
  fit <- update(
    initial_fit,
    newdata = SalmonSurvCUI[1:(i-1), ],
    recompile = FALSE,
    refresh = 0
  )
  
  # Predict the next time point
  new_data <- SalmonSurvCUI[i, , drop = FALSE]
  pred <- posterior_predict(fit, newdata = new_data)
  SalmonSurvCUI$predictions[i - 1] <- mean(pred)
}

plot(SalmonSurvCUI$predictions, SalmonSurvCUI$Survival_scale)
