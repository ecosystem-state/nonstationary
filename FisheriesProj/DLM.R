library(dplyr)
library(MARSS)
library(abind)
juvindices <-readRDS("data/all_juvenile_indices.rds")
con<-readRDS("data/groundfish/wcbts_mean_size_at_age.rds")
con<-readRDS("data/groundfish/wcbts_condition_factor.rds")

unique(con$common_name)
unique(juvindices$common_name)

SalmonSurv <-readRDS("data/Salmon/Survival_covariates.rds")
unique(SalmonSurv$stock_name)
# Pivot the dataframe wider
SalmonSurvCUI<- SalmonSurv%>%filter(stock_name=="Irongate Fall Fingerling")%>%
  dplyr::select(Marine.Survival, brood_year, Beuti_spring, stock_name)#%>%  
 # pivot_wider(names_from = brood_year, values_fill = c(Marine.Survival, Beuti_spring),
  #            values_fill = list(Marine.Survival = NA)) %>%
 # as.data.frame()

## get time indices
years <- SalmonSurvCUI[, 2]
## number of years of data
TT <- length(years)
## get response variable: logit(survival)
dat <- matrix(scale(SalmonSurvCUI[, 1]), nrow = 1)
#get predictor
beuti <- matrix(SalmonSurvCUI[, 3], nrow = 1)
stock <- matrix(as.factor(SalmonSurvCUI[, 4]), nrow = 1)

##get time indices
#m <- dim(beuti)[1] + 1
m<-4

## for process eqn
B <- diag(m)  ## 2x2; Identity
U <- matrix(0, nrow = m, ncol = 1)  ## 2x1; both elements = 0
Q <- matrix(list(0), m, m)  ## 2x2; all 0 for now
#diag(Q) <- c("q.alpha", "q.beta")  ## 2x2; diag = (q1,q2)
diag(Q) <- c("q.alpha", "q.beta", "s.alpha", "s.beta")  ## 2x2; diag = (q1,q2)

## for observation eqn
Z <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z[1, 1, ] <- rep(1, TT)  ## Nx1; 1's for intercept = c
Z[1, 2, ] <- beuti  ## Nx1; predictor variable
Z[1, 3, ] <- 0  ## coefs for stock 2
Z[1, 4, ] <- 0  ## coef for stock 2

Z2 <- array(NA, c(1, m, TT))  ## NxMxT; empty for now
Z2[1, 1, ] <- 0
Z2[1, 2, ] <- 0
Z2[1, 3, ] <- rep(1, TT)  ## Nx1; 1's for intercept = c
Z2[1, 4, ] <- beuti  ## Nx1; predictor variable

Z3<- abind(Z,Z2, along = 1)
#stack arrays of Z on top of each other 
A <- matrix(0)  ## 1x1; scalar = 0
R <- matrix("r")  ## 1x1; scalar = r

## only need starting values for regr parameters
inits_list <- list(x0 = matrix(c(0, 0,0,0), nrow = m))

## list of model matrices & vectors
mod_list <- list(B = B, U = U, Q = Q, Z = Z3, A = A, R = R)

## fit univariate DLM
fit_2 <- MARSS(dat, inits = inits_list, model = mod_list)

## get estimates of alpha
alpha_hat <- c(fit_2$states[1,1],
               fit_2$states[1,-1] - fit_2$states[2,length(dat)])

## get estimates of eta
beta_hat <- fit_2$states[2,]

## plot the estimated level and drift
par(mfrow = c(2,1), mai = c(0.8, 0.8, 0.2, 0.2), omi = c(0, 0, 0, 0))
## plot alpha
plot.ts(alpha_hat, las = 1, lwd = 2, col = "blue",
        ylab = expression(alpha[t]))
## plot eta
plot.ts(beta_hat, las = 1, lwd = 2, col = "blue",
        ylab = expression(beta[t]))
