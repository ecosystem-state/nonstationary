library(dplyr)
library(PNWColors)
library(tidyr)
library(lubridate)
library(readr)
library(stringr)
library(ggplot2)
library(cowplot)
library(Cairo)
library(MCMCvis)
library(HDInterval)
library(reshape2)
library(tidyverse)
library(dplyr)
library(rstan)
library(car)

#####

N <- 90
NP<-3
PN <- rep(seq(1,NP), each = 30)
NR<-3
RN <- rep(seq(1,NR),  30)
K <- 1
x <- matrix(ncol=K, nrow=N)
x[,1]<-rnorm(N, 0,1)
alphaP <- c(2, 4,-3)
betaP<- c(-1, 3, -2)
sigma <-rnorm(N, 1,0.2)
betaR<- c(-0.1, 0.9, -0.8)
y <- c(rep(alphaP, each = N/3)+rep(betaP, each = N/3)*x+rep(betaR, N/3)*x+sigma)


#### STAN data #### 
data <- list(N = N,
             NP=NP,
             P = PN, #population assignment
             NR=NR,
             R = RN, #population assignment
             K = K,# total observations
             x=x,
             y=y
)


#####Assinging Stan Conditions ####

warmups <- 1000
total_iterations <- 3000
max_treedepth <-  10
n_chains <-  3
n_cores <- 4
adapt_delta <- 0.95


####Fitting the Model ####
bhfit <- stan(
  file = here::here("Src/BayesianLinearHierarchicalModels.stan"),
  data = data,
  chains = n_chains,
  warmup = warmups,
  iter = total_iterations,
  cores = n_cores,
  refresh = 250,
  control = list(max_treedepth = max_treedepth,
                 adapt_delta = adapt_delta)
)
MCMCsummary(bhfit)



##### Plotting Results####

simulated.data<- data.frame(y)%>%
  add_column(x = x, period = as.factor(PN), region = as.factor(RN))

col<-pnw_palette("Sunset2",3,type="discrete")

ggplot(data = simulated.data, aes(x = x, y = y, color=period, shape=region)) +
  geom_point(size=0.75,aes(col=period, shape=region)) +
  geom_smooth(method = "lm", se = FALSE, aes(col=period, lty=region)) +
  scale_color_manual(values =  col[1:3])+
  theme_bw()

