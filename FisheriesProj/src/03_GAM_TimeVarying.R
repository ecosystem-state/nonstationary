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
max_tree<-14
adapt_d <- 0.99
iterations<-5000
warmup<-1000
####  No Trends through Time ####
R<-brm(Survival_scale ~ RunType,
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = 3000, cores=3,chains = 3)
conditional_effects(R)
bayes_R2(R)
loo(R)


S<-brm(Survival_scale ~ stock_name,
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(S)
bayes_R2(S)
loo(S)

Reg<-brm(Survival_scale ~ region,
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(Reg)
bayes_R2(Reg)
loo(Reg)

##### Mean Differences; Single trend #####
tR<-brm(Survival_scale ~ RunType+
          s(calendar_year, k=8)+arma(p = 1),
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tR)
bayes_R2(tR)
loo(tR)

tS<-brm(Survival_scale ~ stock_name+
          s(calendar_year, k=8),
         # arma(time=calendar_year, cov = TRUE,gr=stock_name, p = 1),
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=15), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tS)
bayes_R2(tS)
loo(tS)

tReg<-brm(Survival_scale ~ region+
          s(calendar_year, k=8)+
            arma(time=calendar_year,gr=stock_name, p = 1),
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tReg)
bayes_R2(tReg)
loo(tReg)

t<-brm(Survival_scale ~ 
            s(calendar_year, k=8)+
         arma(time=calendar_year,gr=stock_name, p = 1),
          data=SalmonSurvCUI, 
          control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
          iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(t)
bayes_R2(t)
loo(t)

##### Mean Differences; Single trend #####
tRs<-brm(Survival_scale ~
          s(calendar_year,by= RunType, k=8)+
         arma(time=calendar_year,gr=stock_name, p = 1),
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tRs)
bayes_R2(tRs)
loo(tRs)

tSs<-brm(Survival_scale ~ 
          s(calendar_year, by=stock_name, k=8)+
           arma(time=calendar_year,gr=stock_name, p = 1),
        data=SalmonSurvCUI, 
        control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
        iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tSs)
bayes_R2(tSs)
loo(tSs)

tRegs<-brm(Survival_scale ~ 
            s(calendar_year,by = region, k=8),
          data=SalmonSurvCUI, 
          control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
          iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tRegs)
bayes_R2(tRegs)
loo(tRegs)


tC<-brm(Survival_scale ~ stock_name+
             s(calendar_year,by = region, k=8),
           data=SalmonSurvCUI, 
           control=list(adapt_delta=adapt_d,max_treedepth=max_tree), 
           iter = iterations, warmup=warmup,cores=3,chains = 3)
conditional_effects(tC)
bayes_R2(tC)
loo(tC)

