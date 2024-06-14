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

SalmonSurv <- readRDS("data/Salmon/Survival_covariates.rds")
SalmonSurvCUI<- SalmonSurv%>%filter(stock_name=="Columbia River Summers")%>% #grabbing a single stock, I chose at random
  dplyr::select(Marine.Survival, brood_year, Beuti_spring, area,stock_name, RunType)%>% #grabbing columns of interest
  group_by(stock_name)%>%
  mutate(Survival_scale=scale(Marine.Survival), Beuti_scale=scale(Beuti_spring))%>% #making sure data are normalized
  ungroup()%>%
  rename(year=brood_year)%>%
  as.data.frame()

unique(SalmonSurv$stock_name)
unique(SalmonSurv$RunType)
# Try GP model -- this works fine, 3 div transitions 
# fit <- brm(Survival_scale ~ s(year, bs = "gp", k=5) + s(year, bs = "gp", k=5, by = Beuti_scale), data = SalmonSurvCUI,
#                  control=list(adapt_delta=0.99), iter = 3000, chains = 3)
# This gets down to 1 div transition
max_tree<-11
adapt_d <- 0.99


initial_fit <- brm(Survival_scale ~  s(year, k=8) + #s(year,m=2, k=8) +
                     s(Beuti_scale, k=8, by =year),
                    # s(Beuti_spring, area,by =year, k=8, bs="fs", m=2),
                   autocor = cor_ar(~ year, p = 1), 
                   data = SalmonSurvCUI,
                   control=list(adapt_delta=adapt_d,max_treedepth=max_tree), iter = 3000, cores=4,chains = 3)
summ<-summary(initial_fit)
plot(marginal_effects(initial_fit))
bayes_R2(initial_fit)
summ$splines
pp_check(initial_fit)
pp_check(initial_fit, type = "ecdf_overlay")

marginal_effects(initial_fit)
conditions=data.frame(year=c(1990, 2000, 2015))
conditional_effects(initial_fit, conditions = conditions)$Beuti_scale
fit_summary <- summary(initial_fit)
fit_summary$splines

# evaluate fit
pred <- predict(initial_fit)
joined <- cbind(SalmonSurvCUI, pred)
ggplot(joined, aes(year, Estimate)) + 
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5),alpha=0.2) + 
  geom_line() + 
  geom_point(aes(year,Survival_scale)) + 
  ggtitle("Overall fit")


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
