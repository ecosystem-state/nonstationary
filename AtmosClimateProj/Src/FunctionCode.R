##### Writing Bayes DFA model function ####
warmups <- 1000
total_iterations <- 3000
max_treedepth <-  10
n_chains <-  3
n_cores <- 4
adapt_delta <- 0.95
col<-pnw_palette("Sunset2",3,type="discrete")
bayeslinmod <- function(dat, ind,xyz) {
  N <- length(dat$period)
  NP <- 3
  P<- as.numeric(as.factor(dat$period))
  K <- 1
  y <-xyz
  x<-dat[ind]
  data <- list(N = N, #total number of opservations
               NP=NP, #total number of time periods
               P = P, #time period pointer vector
               K = K,#number of covariates, starting with one but can add for final model structure
               x=x, #Upwelling
               y=y #response variable
  )
  
  bhfit <- stan(
    file = here::here("Src/BayesianLinearHierarchicalModels.stan"),
    data = data,
    chains = n_chains,
    warmup = warmups,
    iter = total_iterations,
    cores = n_cores,
    refresh = 250,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta))
  posterior<-data.frame(summary(bhfit, prob=c(0.025, 0.25,0.75, 0.975, 0.1, 0.9))$summary)
  scaled<-NULL
  alphaP<-posterior[1:3,]%>%
    add_column(
      period=rep(c("1983 - 1988","1989 - 2013","2014 - 2022"),1))
  
  betaP<-posterior[4:6,]%>%
    add_column(
      period=rep(c("1983 - 1988","1989 - 2013","2014 - 2022"),1))
  n<- 1000
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled <-rbind(scaled,cbind(tempalpha, tempbeta))
  }
  scaled.anom<<- scaled%>%
    bind_cols(period=rep(rep(c(1,2, 3), each = 1000),1), Index=rep(colnames(dat[ind]), 3000),
              yfirst=rep(min(dat$Year_lag), 3000),ylast=rep(max(dat$Year_lag), 3000))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}

postplot <- function(dat, param){
  ggplot(dat, aes(x = param, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1], col[2], col[3])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
}

postplot(CALCOFI, CALCOFI$beta)
##### DFA Trend Model Runs ####
season<- "Winter"
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%rename(period=era)%>%
  mutate(period=ifelse(Year_lag<=1988,1, ifelse(Year_lag>2004,3,2)))
dat <- dfa%>%select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

data <- dat%>%filter(season==season, trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)

CALCOFI <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1)
  CALCOFI <-rbind(CALCOFI,scaled.anom)
}
CALCOFI<-CALCOFI%>%mutate(survey="CALCOFI")
postplot(CALCOFI, CALCOFI$beta)
postplot(CALCOFI, CALCOFI$alpha)


data <- dat%>%filter(season==season, trend=="SEA")%>%
  distinct()%>%
  filter(Year_lag<2023)
SEA <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1)
  SEA <-rbind(SEA,scaled.anom)
}
SEA<-SEA%>%mutate(survey="SEA")
postplot(SEA, SEA$beta)
postplot(SEA, SEA$alpha)


data <- dat%>%filter(season==season, trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)

RREAS <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1)
  RREAS <-rbind(RREAS,scaled.anom)
}
RREAS<-RREAS%>%mutate(survey="RREAS")
postplot(RREAS, RREAS$beta)
postplot(RREAS, RREAS$alpha)

#full.dfa.spring<-RREAS%>%bind_rows(SEA,CALCOFI)

full.dfa.winter.offset<-RREAS%>%bind_rows(SEA,CALCOFI)

dfa.spring.slope<- ggplot(full.dfa.winter.offset, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
    facet_grid(Index~survey, scales='free') +
    geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Figures/survey.winter.slope.offset.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8)
dfa.spring.slope
dev.off()

#### Upwelling Model Runs ####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
data<- climate_dat%>%filter(season=="Winter")%>%
    distinct()%>%
    filter(Year_lag<2023)

bayeslinmod_up <- function(dat, ind,xyz) {
  N <- length(dat$period)
  NP <- 12
  P<- as.numeric(as.factor(dat$era.region))
  K <- 1
  y <-xyz
  x<-dat[ind]
  data <- list(N = N, #total number of opservations
               NP=NP, #total number of time periods
               P = P, #time period pointer vector
               K = K,#number of covariates, starting with one but can add for final model structure
               x=x, #Upwelling
               y=y #response variable
  )
  warmups <- 1000
  total_iterations <- 3000
  max_treedepth <-  10
  n_chains <-  3
  n_cores <- 4
  adapt_delta <- 0.95
  
  bhfit <- stan(
    file = here::here("Src/BayesianLinearHierarchicalModels.stan"),
    data = data,
    chains = n_chains,
    warmup = warmups,
    iter = total_iterations,
    cores = n_cores,
    refresh = 250,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta))
  
  posterior<-data.frame(summary(bhfit, prob=c(0.025, 0.25,0.75, 0.975, 0.1, 0.9))$summary)
  alphaP<-posterior[1:12,]%>%
    add_column(region = rep(c("GoA", "NCC", "CCC","SCC"),each =3), 
               period=rep(c("1967 - 1988","1989 - 2013","2014 - 2022"),4))
  
  betaP<-posterior[13:24,]%>%
    add_column(region = rep(c("GoA", "NCC","CCC", "SCC"),each =3), 
               period=rep(c("1967 - 1988","1989 - 2013","2014 - 2022"),4))
  n<- 1000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  scaled<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('1967 - 1988', '1989 - 2013', '2014 - 2022'), each = 1000),4),
              region = rep(rep(c("GoA", "NCC","CCC", "SCC"),each = 3000)), Index = rep(colnames(dat[ind]), 12000))%>%
    rename(alpha=tempalpha, beta = tempbeta)

}

columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
up_post<-NULL
for(i in 1:length(columns)){
  bayeslinmod_up(data, columns[i], data$stand_bakun_seasonally)
  up_post <-rbind(up_post,scaled)
}

region.lev2 <-c("GoA", "NCC", "CCC", "SCC")
up_post$region <-factor(up_post$region, levels=region.lev2)
up_post_winter<-up_post
up.spring.int<-ggplot(up_post_spring, aes(x = alpha, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

up.spring.slope<-ggplot(up_post_spring, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Figures/up.spring.slope.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8)
up.spring.slope
dev.off()

#### Copepod Model Runs ####
bayeslinmod_cop <- function(dat, ind,xyz) {
  N <- length(dat$period)
  NP <- 2
  P<- as.numeric(as.factor(dat$period))
  K <- 1
  y <-xyz
  x<-dat[ind]
  data <- list(N = N, #total number of opservations
               NP=NP, #total number of time periods
               P = P, #time period pointer vector
               K = K,#number of covariates, starting with one but can add for final model structure
               x=x, #Upwelling
               y=y #response variable
  )
  
  bhfit <- stan(
    file = here::here("Src/BayesianLinearHierarchicalModels.stan"),
    data = data,
    chains = n_chains,
    warmup = warmups,
    iter = total_iterations,
    cores = n_cores,
    refresh = 250,
    control = list(max_treedepth = max_treedepth,
                   adapt_delta = adapt_delta))
  posterior<<-data.frame(summary(bhfit, prob=c(0.025, 0.25,0.75, 0.975, 0.1, 0.9))$summary)
  
  alphaP<-posterior[1:2,]%>%
    add_column(
      period=rep(c("1989 - 2013","2014 - 2022"),1))
  
  betaP<-posterior[3:4,]%>%
    add_column(
      period=rep(c("1989 - 2013","2014 - 2022"),1))
  n<- 1000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  
  scaled.anom<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('1989 - 2013', '2014 - 2022'), each = 1000),1), Index = rep(colnames(dat[ind]), 2000))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}
climpdo <- climate_dat_cop%>%filter(season=="Winter")%>%
  select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Spring")%>%
  select(Year_lag, region, seasonal_copepod_northern, annual_copepod_northern, 
         seasonal_copepod_southern)%>%
  filter(Year_lag<2023&Year_lag>1996)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))

northern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_northern)
  northern <-rbind(northern,scaled.anom)
}
northern<-northern%>%mutate(survey="northern")

southern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_southern)
  southern <-rbind(southern,scaled.anom)
}
southern<-southern%>%mutate(survey="southern")
COP<-southern%>%bind_rows(northern)


cop.wisp.slope<-ggplot(COP, aes(x = beta, fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[2], col[3])) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Figures/cop.wisp.slope.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8)
cop.wisp.slope
dev.off()

