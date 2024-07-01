### Model 1: time varying intercept and slope
data(SalmonSurvCUI)
head(SalmonSurvCUI)

SalmonSurvCUI$CUI.apr = scale(SalmonSurvCUI$CUI.apr)

fit <- fit_dlm(time_varying = logit.s ~ CUI.apr,
        data = SalmonSurvCUI,
        chains=1,
        iter=1000)

fit <- fit_dlm(time_varying = logit.s ~ CUI.apr,
        data = SalmonSurvCUI[-36,],
        chains=1,
        iter=1000)



fit <- fit_dlm(time_varying = logit.s ~ CUI.apr,
        data = SalmonSurvCUI,
        chains=1,
        iter=1000,
        control = list(adapt_delta=0.99))



fit <- fit_dlm(time_varying = logit.s ~ CUI.apr,
        data = SalmonSurvCUI,
        chains=3,
        iter=12000,
        warmup = 3000,
        control = list(max_treedepth=15, adapt_delta=0.99))



pairs(fit)
dlm_trends(fit)



### Model 2: time varying intercept and constant slope

fit <- fit_dlm(time_varying = logit.s ~ 1,
               formula = logit.s ~ CUI.apr,
        data = SalmonSurvCUI,
        chains=3,
        iter=10000,control = list(max_treedepth=15, adapt_delta=0.99),
        warmup = 3000)


### Model 3: constant intercept and time varying slope

fit <- fit_dlm(time_varying = logit.s ~ CUI.apr,
               formula = logit.s ~ 1,
        data = SalmonSurvCUI,
        chains=1,
        iter=10000,
        warmup = 3000,
        control = list(max_treedepth=15, adapt_delta=0.99))


### Comparing models


library(loo)
loo::loo(fit$fit)