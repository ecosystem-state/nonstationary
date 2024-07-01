library(dplyr)
library(MARSS)
library(abind)
library(RcppEigen)
library(StanHeaders)
library(mvtnorm)
#remotes::install_github("atsa-es/mvdlm")
library(mvdlm)
library(ggplot2)
library(broom.mixed)
library(loo)

#### Reading in Data ####
SalmonSurv <-readRDS("data/Salmon/Survival_covariates.rds") ###reading in salmon data
unique(SalmonSurv$stock_name) #looking up stock names to pull one out

SalmonSurvCUI<- SalmonSurv%>%filter(stock_name=="Atnarko")%>% #grabbing a single stock, I chose at random
  dplyr::select(Marine.Survival, brood_year, Beuti_spring, stock_name)%>% #grabbing columns of interest
  mutate(Survival_scale=scale(Marine.Survival), Beuti_scale=scale(Beuti_spring))%>% #making sure data are normalized
  rename(year=brood_year)%>%
 as.data.frame()

#### STAN controls ####
#chains <- 3
chains<-1
iterations<-10000
warmups <- 3000
max_tree<-15
adapt_d <- 0.99

#### Time-varying intercept, constant slope ####

fit_intercept <- fit_dlm(time_varying = Survival_scale  ~ 1,
               formula = Survival_scale ~ Beuti_scale,
        data = SalmonSurvCUI,
        chains=chains,
        iter=iterations,
        warmup = warmups,
        control = list(max_treedepth=max_tree, adapt_delta=adapt_d))

#Warnings: 105 divergent transitions; BFMI is low
dlm_trends(fit_intercept) 

#### Time-varying slope; constant intercept ####

fit_slope <- fit_dlm(time_varying = Survival_scale ~ Beuti_scale,
               formula = Survival_scale ~ 1,
               data = SalmonSurvCUI,
               chains=chains,
               iter=iterations,
               warmup = warmups,
               control = list(max_treedepth=max_tree, adapt_delta=adapt_d))

#Warnings: 2271 divergent transitions; BFMI low, R-hat is high and ESS is low
dlm_trends(fit_slope)#when I plot this it still looks like slope and intercept are both time-varying????

#### Time-varying intercept and time-varying slope ####

fit_both <- fit_dlm(time_varying = Survival_scale  ~ Beuti_scale,
               data = SalmonSurvCUI,
               chains=chains,
               iter=iterations,
               warmup = warmups,
               correlated_rw=TRUE,
               control = list(max_treedepth=max_tree, adapt_delta=adapt_d))

broom.mixed::tidy(fit_both$fit)
dlm_trends(fit_both)

#Warnings: 25 divergent transitions, BFMI is low...not TOO bad...

#### Running with correlated RW ####

fit_both <- fit_dlm(time_varying = Survival_scale  ~ Beuti_scale,
                    data = SalmonSurvCUI,
                    chains=chains,
                    iter=iterations,
                    warmup = warmups,
                    correlated_rw=TRUE,
                    control = list(max_treedepth=max_tree, adapt_delta=adapt_d))

broom.mixed::tidy(fit_both$fit)
dlm_trends(fit_both)

#Warnings 525 divergent transtins, BFMI is low; ESS is low
