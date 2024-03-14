library(ncdf4)
library(chron)
library(tidyverse)
library(kohonen) # fitting
library(aweSOM) # plotting
library(SOMbrero) # plotting
library(paletteer) #colors
library(PNWColors) #more colors
library(here) #navigating folders
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(maps)       #basic mapping functions and some data
library(rstan)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(PBSmapping)
library(overlapping)
set.seed(1234)

##### Writing Bayes DFA model function ####
warmups <- 1000
total_iterations <- 10000
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
  n<- 10000
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled <-rbind(scaled,cbind(tempalpha, tempbeta))
  }
  scaled.anom<<- scaled%>%
    bind_cols(period=rep(rep(c(1,2, 3), each = n),1), Index=rep(colnames(dat[ind]), n*3),
              yfirst=rep(min(dat$Year_lag), n*3),ylast=rep(max(dat$Year_lag), n*3))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}

postplot <- function(dat, param){
  ggplot(dat, aes(x = param, fill = as.factor(period), group=as.factor(period))) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1], col[2], col[3], "grey")) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")
}

##### DFA Trend Model Runs ####
climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

climpdo <- climate_dat_cop%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)
data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
         seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  pivot_longer(!c(Year_lag,  NPH, NPGO,PDO,ONI), 
               names_to = "trend", values_to = "estimate")%>%
  mutate(season="Spring")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.long<- dat%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, estimate,trend,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.lm<-data%>%
  bind_rows(dat.long)%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  filter(trend!="SEA")%>%
  mutate(trend=ifelse(trend=="CALCOFI","CALCOFI (SCC)",
                      ifelse(trend=="RREAS", "RREAS (CCC)", 
                          ifelse(trend=="seasonal_copepod_southern","Southern Copepod (NCC)","Northern Copepod (NCC)"))))%>%
  mutate(trend = fct_relevel(trend, "Northern Copepod (NCC)", 
                              "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))

survey.lm<-ggplot(data = dat.lm, aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
   geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")

pdf(file = "Output/Figures/biological.lm.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
survey.lm
dev.off()
###### CALCOFI ######

data <- dat%>%filter(season=="Spring", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)

CALCOFI <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate)
  CALCOFI <-rbind(CALCOFI,scaled.anom)
}
CALCOFI<-CALCOFI%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI, CALCOFI$beta)
postplot(CALCOFI, CALCOFI$alpha)

ratio.beta <- cbind(CALCOFI.beta=c(CALCOFI%>%filter(period=='3')%>%dplyr::select(beta)/CALCOFI%>%filter(period=='2')%>%dplyr::select(beta)),
                    CALCOFI%>%filter(period=='3')%>%dplyr::select(Index))%>%
  rename(CALCOFI=beta)
ratio.alpha<-cbind(CALCOFI.alpha=CALCOFI%>%filter(period=='3')%>%dplyr::select(alpha)/CALCOFI%>%filter(period=='2')%>%dplyr::select(alpha), 
                   CALCOFI%>%filter(period=='3')%>%dplyr::select(Index))%>%
  rename(CALCOFI=alpha)

index.names <- unique(CALCOFI$Index)
period.names <- unique(CALCOFI$period)
overlap.CALCOFI <- NA
  for(i in 1:4){
    temp <- CALCOFI%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.CALCOFI <-rbind(temp2,overlap.CALCOFI)
  }
overlap.CALCOFI 

#running models with 1-year offset

CALCOFIoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1)
  CALCOFIoffset <-rbind(CALCOFIoffset,scaled.anom)
}
CALCOFIoffset<-CALCOFIoffset%>%mutate(survey="CALCOFI (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(CALCOFIoffset, CALCOFIoffset$beta)
postplot(CALCOFIoffset, CALCOFIoffset$alpha)


###### RREAS ######

data <- dat%>%filter(season=="Spring", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)

RREAS <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate)
  RREAS <-rbind(RREAS,scaled.anom)
}
RREAS<-RREAS%>%mutate(survey="RREAS (CCC)",
                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                          ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  filter(period!=1)
  
postplot(RREAS, RREAS$beta)
postplot(RREAS, RREAS$alpha)

index.names <- unique(RREAS$Index)
period.names <- unique(RREAS$period)
overlap.RREAS <- NA
for(i in 1:4){
  temp <- RREAS%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.RREAS <-rbind(temp2,overlap.RREAS)
}
overlap.RREAS
overlap.CALCOFI

RREASoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1)
  RREASoffset <-rbind(RREASoffset,scaled.anom)
}
RREASoffset<-RREASoffset%>%mutate(survey="RREAS (CCC)",
                                  Index=ifelse(Index=="seasonal_PDO", "PDO", 
                                      ifelse(Index=="seasonal_NPGO", "NPGO",
                                      ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  filter(period!=1)
postplot(RREASoffset, RREASoffset$beta)
postplot(RREASoffset, RREASoffset$alpha)

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
      period=rep(c("2","3"),1))
  
  betaP<-posterior[3:4,]%>%
    add_column(
      period=rep(c("2","3"),1))
  n<- 10000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  
  scaled.anom<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('2', '3'), each = n),1), Index = rep(colnames(dat[ind]), n*2))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}


climpdo <- climate_dat_cop%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
         seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)

northern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_northern)
  northern <-rbind(northern,scaled.anom)
}
northern<-northern%>%mutate(survey="Northern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))



index.names <- unique(northern$Index)
period.names <- unique(northern$period)


overlap.northern <- NA
for(i in 1:4){
  temp <- northern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.northern <-rbind(temp2,overlap.northern)
}
overlap.northern

southern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_southern)
  southern <-rbind(southern,scaled.anom)
}
southern<-southern%>%mutate(survey="Southern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
index.names <- unique(southern$Index)
period.names <- unique(southern$period)
overlap.southern <- NA
for(i in 1:4){
  temp <- southern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.southern <-rbind(temp2,overlap.southern)
}
overlap.southern


COP<-southern%>%bind_rows(northern)%>%mutate(yfirst=1989, ylast=2023,period=as.numeric(period))


####### Combining data for full plots #####
#full.dfa.spring<-RREAS%>%bind_rows(SEA,CALCOFI)

fulldat<-RREAS%>%bind_rows(CALCOFI, COP)%>%
  mutate(survey = fct_relevel(survey, "Northern Copepod (NCC)", 
                              "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"),
         period2 = ifelse(period==1,"1967 - 1988", ifelse(period==2, "1989 - 2012", "2013 - 2022")))

full.datslope<- ggplot(fulldat, aes(x = beta, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
    facet_grid(Index~survey, scales='free') +
    geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Figures/biological.slope.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
full.datslope
dev.off()

full.datint<- ggplot(fulldat, aes(x = alpha, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

pdf(file = "Output/Figures/biological.intercept.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
full.datint
dev.off()


fulloffset<-RREASoffset%>%bind_rows(CALCOFIoffset)%>%
  mutate(survey = fct_relevel(survey, "RREAS (CCC)","CALCOFI (SCC)"),
         period2 = ifelse(period==1,"1967 - 1988", ifelse(period==2, "1989 - 2012", "2013 - 2022")))


offset.slope<- ggplot(fulloffset, aes(x = beta, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-Year Lag")

pdf(file = "Output/Figures/biologicaloffset.slope.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
offset.slope
dev.off()

offset.int<- ggplot(fulloffset, aes(x = alpha, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-Year Lag")
pdf(file = "Output/Figures/biologicaloffset.intercept.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
offset.int
dev.off()


#### Upwelling Model Runs ####

##### spring #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Spring"
data<- climate_dat%>%filter(season=="Spring")%>%
#  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()
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
  total_iterations <- 10000
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
  n<- 10000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  scaled<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('1967 - 1988', '1989 - 2013', '2014 - 2022'), each = n),4),
              region = rep(rep(c("GoA", "NCC","CCC", "SCC"),each = n*3)), Index = rep(colnames(dat[ind]), 12*n))%>%
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
up_post_spring<-up_post%>%
  mutate(era=as.numeric(as.factor(period)),
         period2 = ifelse(era==1,"1967 - 1988", ifelse(era==2, "1989 - 2012", "2013 - 2023")),
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
index.names <- unique(up_post_spring$Index)
region.names<-unique(up_post_spring$region)
period.names <- unique(up_post_spring$period)
overlap.up <- NA
for(j in 1:4){
  temp1 <- up_post_spring%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period=='1967 - 1988')%>%dplyr::select(beta),temp%>%filter(period=="1989 - 2013")%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period=='1967 - 1988')%>%dplyr::select(beta),temp%>%filter(period=="2014 - 2022")%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period=="2014 - 2022")%>%dplyr::select(beta),temp%>%filter(period=="1989 - 2013")%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.up <-rbind(temp2,overlap.up)
  }

}

overlap.up 
mean(na.omit(overlap.up%>%filter(period1==1&period2==3)%>%dplyr::select(ov))$ov)
sd(na.omit(overlap.up%>%filter(period1==1&period2==3)%>%dplyr::select(ov))$ov)

mean(na.omit(overlap.up%>%filter(period1==3&period2==2)%>%dplyr::select(ov))$ov)
sd(na.omit(overlap.up%>%filter(period1==3&period2==2)%>%dplyr::select(ov))$ov)

up.spring.int<-ggplot(up_post_spring%>%filter(region!="GoA"), aes(x = alpha, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")+
  ggtitle("Spring")

up.spring.slope<-ggplot(up_post_spring%>%filter(region!="GoA"), aes(x = beta, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")+
  ggtitle("Spring")

pdf(file = "Output/Figures/up.slope.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
up.spring.slope
dev.off()

pdf(file = "Output/Figures/up.int.pdf",   # The directory you want to save the file in
    width =8, # The width of the plot in inches
    height = 8)
up.spring.int
dev.off()


dat.lm<- data%>%dplyr::select(Year_lag, period,region, season, era.region,stand_bakun_seasonally,seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, period, region, era.region,season, stand_bakun_seasonally), names_to = "Index_Name", values_to = "Index_Value")
region.lev3 <-c("GoA", "Northern CC", "Central CC", "Southern CC")
dat.lm$region <-factor(dat.lm$region, levels=region.lev3)

up.lm<-ggplot(data = dat.lm%>%filter(region!="GoA"), aes(y = stand_bakun_seasonally, x =Index_Value,col=period)) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=period)) +
  #  geom_text(aes(label=Year_lag,col=period)) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  scale_y_continuous(name = "Upwelling (Bakun 1º 6-hourly)") +
  scale_color_manual(values =  col[1:3], labels=c('1967 - 1988', '1989 - 2012', '2013 - 2023'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")

pdf(file = "Output/Figures/up.lm.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8.5)
up.lm
dev.off()

##### winter#####

climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Winter"
data<- climate_dat%>%filter(season=="Winter")%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()
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
  total_iterations <- 10000
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
  n<- 10000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  scaled<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('1967 - 1988', '1989 - 2013', '2014 - 2022'), each = n),4),
              region = rep(rep(c("GoA", "NCC","CCC", "SCC"),each = n*3)), Index = rep(colnames(dat[ind]), 12*n))%>%
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
up_post_winter<-up_post%>%
  mutate(era=as.numeric(as.factor(period)),
         period2 = ifelse(era==1,"1967 - 1988", ifelse(era==2, "1989 - 2012", "2013 - 2023")),
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
index.names <- unique(up_post_winter$Index)
region.names<-unique(up_post_winter$region)
period.names <- unique(up_post_winter$period)
overlap.up <- NA
for(j in 1:4){
  temp1 <- up_post_winter%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period=='1967 - 1988')%>%dplyr::select(beta),temp%>%filter(period=="1989 - 2013")%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period=='1967 - 1988')%>%dplyr::select(beta),temp%>%filter(period=="2014 - 2022")%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period=="2014 - 2022")%>%dplyr::select(beta),temp%>%filter(period=="1989 - 2013")%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    temp2<-rbind(ov1,ov2,ov3)
    overlap.up <-rbind(temp2,overlap.up)
  }
  
}

overlap.up 
mean(na.omit(overlap.up%>%filter(period1==1&period2==3)%>%dplyr::select(ov))$ov)
sd(na.omit(overlap.up%>%filter(period1==1&period2==3)%>%dplyr::select(ov))$ov)

mean(na.omit(overlap.up%>%filter(period1==3&period2==2)%>%dplyr::select(ov))$ov)
sd(na.omit(overlap.up%>%filter(period1==3&period2==2)%>%dplyr::select(ov))$ov)

up.spring.int<-ggplot(up_post_winter%>%filter(region!="GoA"), aes(x = alpha, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")+
  ggtitle("Winter")

up.spring.slope<-ggplot(up_post_winter%>%filter(region!="GoA"), aes(x = beta, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~region, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  theme(plot.title = element_text(hjust = 0.5))+
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")+
  ggtitle("Winter")

pdf(file = "Output/Figures/up.slope.winter.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
up.spring.slope
dev.off()

pdf(file = "Output/Figures/up.int.winter.pdf",   # The directory you want to save the file in
    width =8, # The width of the plot in inches
    height = 8)
up.spring.int
dev.off()


dat.lm<- data%>%dplyr::select(Year_lag, period,region, season, era.region,stand_bakun_seasonally,seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, period, region, era.region,season, stand_bakun_seasonally), names_to = "Index_Name", values_to = "Index_Value")
region.lev3 <-c("GoA", "Northern CC", "Central CC", "Southern CC")
dat.lm$region <-factor(dat.lm$region, levels=region.lev3)

up.lm<-ggplot(data = dat.lm%>%filter(region!="GoA"), aes(y = stand_bakun_seasonally, x =Index_Value,col=period)) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=period)) +
  #  geom_text(aes(label=Year_lag,col=period)) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  scale_y_continuous(name = "Upwelling (Bakun 1º 6-hourly)") +
  scale_color_manual(values =  col[1:3], labels=c('1967 - 1988', '1989 - 2012', '2013 - 2023'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")

pdf(file = "Output/Figures/up.lm.winter.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8.5)
up.lm
dev.off()
#### Violin plot Uowelling data file #####

ratio.beta.up<-rbind(up_post%>%filter(period=='1967 - 1988')%>%add_column((up_post%>%filter(period=='2014 - 2022')%>%group_by(Index,region))$beta-
                                                                         (up_post%>%filter(period=='1989 - 2013')%>%group_by(Index,region))$beta)%>%dplyr::select(region,beta,Index,`... - ...`)%>%
                                                                         rename('beta_diff'=`... - ...`)%>%
                       mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2"),
                  up_post%>%filter(period=='1967 - 1988')%>%add_column((up_post%>%filter(period=='2014 - 2022')%>%group_by(Index,region))$beta-
                                                                         (up_post%>%filter(period=='1967 - 1988')%>%group_by(Index,region))$beta)%>%dplyr::select(region,beta,Index,`... - ...`)%>%
                    rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - 1967:1988",Difference2="Era 3 - Era 1"),
                  up_post%>%filter(period=='1967 - 1988')%>%add_column((up_post%>%filter(period=='1989 - 2013')%>%group_by(Index,region))$beta-
                                                                         (up_post%>%filter(period=='1967 - 1988')%>%group_by(Index,region))$beta)%>%dplyr::select(region,beta,Index,`... - ...`)%>%
                    rename('beta_diff'=`... - ...`)%>%mutate(Difference="1989:2012 - 1967:1988", Difference2="Era 2 - Era 1"))%>%
  mutate(region = fct_relevel(region, "SCC","CCC", "NCC"))


# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }

col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ratio.up<-ratio.beta.up%>%filter(region!="GoA")%>%
mutate(Index2=ifelse(Index=='seasonal_NPGO','NPGO', ifelse(Index=='seasonal_NPH','NPH', ifelse(Index=='seasonal_ONI', 'ONI','PDO'))))

violin.up <-ggplot(ratio.up, aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index2~Difference2) +
  ylab("Upwelling Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.up 
pdf(file = "Output/Figures/violin.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
violin.up
dev.off()


#### Violin plot BIOLOGICAL data file #####

ratio.beta<-rbind(CALCOFI%>%filter(period=='3')%>%add_column((CALCOFI%>%filter(period=='3')%>%group_by(Index))$beta-(CALCOFI%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
northern%>%filter(period=='3')%>%add_column((northern%>%filter(period=='3')%>%group_by(Index))$beta-(northern%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
southern%>%filter(period=='3')%>%add_column((southern%>%filter(period=='3')%>%group_by(Index))$beta-(southern%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
RREAS%>%filter(period=='3')%>%add_column((RREAS%>%filter(period=='3')%>%group_by(Index))$beta-(RREAS%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`))

colnames(ratio.beta)<- c("alpha", "beta", "survey","Index", "ratio")

ratio.alpha<-rbind(CALCOFI%>%filter(period=='3')%>%add_column((CALCOFI%>%filter(period=='3')%>%group_by(Index))$alpha-(CALCOFI%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  northern%>%filter(period=='3')%>%add_column((northern%>%filter(period=='3')%>%group_by(Index))$alpha-(northern%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  southern%>%filter(period=='3')%>%add_column((southern%>%filter(period=='3')%>%group_by(Index))$alpha-(southern%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  RREAS%>%filter(period=='3')%>%add_column((RREAS%>%filter(period=='3')%>%group_by(Index))$alpha-(RREAS%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`))

colnames(ratio.alpha)<- c("alpha", "beta", "survey", "Index","ratio")

# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }
ratio.alpha<-ratio.alpha%>%
  mutate(parameter="Intercept")

ratio.beta<-ratio.beta%>%
  mutate(parameter="Slope")

ratio<-ratio.beta%>%#add_row(ratio.alpha)%>%
  mutate(Difference="2013:2023 - 1989:2012",Difference2="Era 3 - Era 2")%>%
  mutate(survey2 = ifelse(survey=="CALCOFI (SCC)", "CalCOFI",
                          ifelse(survey=="RREAS (CCC)", "RREAS",
                                 ifelse(survey=="Southern Copepod (NCC)","s. Copepod","n. Copepod"))))%>%
  mutate(survey2 = factor(survey2, levels=c(c("CalCOFI","RREAS", "s. Copepod","n. Copepod"))))
  
col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

violin <-ggplot(ratio, aes(x=survey2, y=ratio, fill=survey2)) +
  theme_bw() +
  scale_fill_manual(values=c(col4[3:1],col4[1]), name="Survey")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Difference2~Index) +
  ylab("Biological Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) + theme(legend.position = "none")
violin
pdf(file = "Output/Figures/violin.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
violin
dev.off()

#### Map Data #####
col2<-pnw_palette("Sunset2",4,type="discrete")
col<-pnw_palette("Sunset2",3,type="discrete")
col3<-pnw_palette("Sunset2",8,type="continuous")

col<-pnw_palette("Sunset2",3,type="discrete")
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
climate_dat_cop <-readRDS(here('data/physical/climate_dat_cop.rds'))
bakunsites <- read.csv(here('data/physical/Bakun/MapLocations.csv'))%>%
  mutate(longitude=longitude)
sites <- st_as_sf(data.frame(bakunsites[,1:2]), coords = c("longitude","latitude"), crs = 4326, 
                  agr = "constant")

#### Making the Map #####
map<-ggplot() +
  geom_polygon(aes(x=c(-105, -113, -127,-105,-105),
                   y=c(22.1, 22.1,34.4486,34.4486,20)),
               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110, -127,-130,-110,-110), 
                   y=c(34.4486,34.4486,40.4401,40.4401,34.4486)),
               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110, -130,-130,-110,-110), 
                   y=c(40.4401,40.4401,49.5,49.5,40.4401)),
               fill='white', col='black',alpha=0.6)+
  geom_sf(data = world)+
  annotate("rect", xmin= -121.5, xmax = -109, ymin = 42, ymax = 48.8, 
           fill = 'white', col='black',size = 0.8, lwd=0.2) +
  geom_sf(fill='grey95') +
  geom_sf(data = sites, size = c(rep(2,68+35+16), rep(3,2)), 
          shape = c(rep(24,68), rep(21,35),rep(23,16),rep(22,2)), 
          col = c(rep('black',68+35+18)), 
          fill = c(rep(col3[3],68), rep(col3[7],35),rep(col3[5],16),rep(col3[8],2))) +
  coord_sf(xlim = c(-132, -108), ylim = c(22, 50), expand = FALSE)+
  ylab(" ")+
  xlab(" ")+
  annotation_scale()+
  annotation_north_arrow(which_north = "true",pad_x = unit(0.25, "in"), 
                         pad_y = unit(0.25, "in"))+
  annotate(geom = "text", x = c(-127,-118.5,-129), y = c(37,28,45), 
           label = str_wrap(c("Central", "Southern","Northern"), width = 20),
           fontface = "italic", color = "grey22", size = 3.75, angle=c('290', '311','270')) +
  annotate(geom = "text", x = c(-114,-114,-114,-114,-114), y = c(48,46.5,45, 44,43), 
           label = str_wrap(c("Bakun Index","CC Regions", "Newport Line","CalCOFI", "RREAS"), width = 22),
           color = "grey22", size =3.5) +
  annotate(geom = "text", x = c(-117,-113), y = c(41,35), 
           label = str_wrap(c("Cape Mendocino","Point Conception"), width = 20),
           fontface = "italic", color = "grey22", size = 3) +
  annotate("rect", xmin= -121, xmax = -119, ymin = 46, ymax = 47, 
           fill = 'white', col='black',size = 0.8, lwd=0.5) +
  annotate("line", x= c(-124.1, -124.65), y = c(44.652, 44.652),col=col2[1],size = 0.8, lwd=1) +
  annotate("line", x= c(-120.5, -119.5), y = c(45, 45),col=col2[1],size = 0.8, lwd=1) +
  
  theme(panel.background = element_rect(fill = "lightsteelblue2"),
        panel.border = element_rect(fill = NA),panel.grid.major = element_line(colour = "transparent"))

map 
#### JUMBO Plot #####
pdf("Output/Fig Jumbo.pdf", 11,8) 
ggarrange(ggarrange(map, violin.up, ncol = 2, labels = c("A", "B")),violin,  # Second row with box and dot plots
          nrow = 2, heights=c(2,0.75),
          labels = c("A", "C")                                        # Labels of the scatter plot
) 
dev.off()



##### Writing Bayes DFA model function ####

##### DFA Trend Model Runs ####
climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

climpdo <- climate_dat_cop%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)
data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
         seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  pivot_longer(!c(Year_lag,  NPH, NPGO,PDO,ONI), 
               names_to = "trend", values_to = "estimate")%>%
  mutate(season="Spring")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.long<- dat%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, estimate,trend,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.lm<-data%>%
  bind_rows(dat.long)%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  filter(trend!="SEA")%>%
  mutate(trend=ifelse(trend=="CALCOFI","CALCOFI (SCC)",
                      ifelse(trend=="RREAS", "RREAS (CCC)", 
                             ifelse(trend=="seasonal_copepod_southern","Southern Copepod (NCC)","Northern Copepod (NCC)"))))%>%
  mutate(trend = fct_relevel(trend, "Northern Copepod (NCC)", 
                             "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))

survey.lm<-ggplot(data = dat.lm, aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle(season)

pdf(file = "Output/Figures/biological.winter.lm.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
survey.lm
dev.off()
###### CALCOFI ######

data <- dat%>%filter(season=="Winter", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)

CALCOFI <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate)
  CALCOFI <-rbind(CALCOFI,scaled.anom)
}
CALCOFI<-CALCOFI%>%mutate(survey="CALCOFI (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI, CALCOFI$beta)
postplot(CALCOFI, CALCOFI$alpha)

ratio.beta <- cbind(CALCOFI.beta=c(CALCOFI%>%filter(period=='3')%>%dplyr::select(beta)/CALCOFI%>%filter(period=='2')%>%dplyr::select(beta)),
                    CALCOFI%>%filter(period=='3')%>%dplyr::select(Index))%>%
  rename(CALCOFI=beta)
ratio.alpha<-cbind(CALCOFI.alpha=CALCOFI%>%filter(period=='3')%>%dplyr::select(alpha)/CALCOFI%>%filter(period=='2')%>%dplyr::select(alpha), 
                   CALCOFI%>%filter(period=='3')%>%dplyr::select(Index))%>%
  rename(CALCOFI=alpha)

index.names <- unique(CALCOFI$Index)
period.names <- unique(CALCOFI$period)
overlap.CALCOFI <- NA
for(i in 1:4){
  temp <- CALCOFI%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
  ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
  temp2<-rbind(ov1,ov2,ov3)
  overlap.CALCOFI <-rbind(temp2,overlap.CALCOFI)
}
overlap.CALCOFI 

###### RREAS ######

data <- dat%>%filter(season=="Winter", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)

RREAS <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate)
  RREAS <-rbind(RREAS,scaled.anom)
}
RREAS<-RREAS%>%mutate(survey="RREAS (CCC)",
                      Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                        ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  filter(period!=1)

postplot(RREAS, RREAS$beta)
postplot(RREAS, RREAS$alpha)

index.names <- unique(RREAS$Index)
period.names <- unique(RREAS$period)
overlap.RREAS <- NA
for(i in 1:4){
  temp <- RREAS%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.RREAS <-rbind(temp2,overlap.RREAS)
}
overlap.RREAS
overlap.CALCOFI



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
      period=rep(c("2","3"),1))
  
  betaP<-posterior[3:4,]%>%
    add_column(
      period=rep(c("2","3"),1))
  n<- 10000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  
  scaled.anom<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('2', '3'), each = n),1), Index = rep(colnames(dat[ind]), n*2))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}


climpdo <- climate_dat_cop%>%filter(season=="Winter")%>%
  dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
         seasonal_copepod_southern)%>%
  left_join(climpdo)%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)

northern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_northern)
  northern <-rbind(northern,scaled.anom)
}
northern<-northern%>%mutate(survey="Northern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))



index.names <- unique(northern$Index)
period.names <- unique(northern$period)


overlap.northern <- NA
for(i in 1:4){
  temp <- northern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.northern <-rbind(temp2,overlap.northern)
}
overlap.northern

southern <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod_cop(data, columns[i],  data$seasonal_copepod_southern)
  southern <-rbind(southern,scaled.anom)
}
southern<-southern%>%mutate(survey="Southern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
index.names <- unique(southern$Index)
period.names <- unique(southern$period)
overlap.southern <- NA
for(i in 1:4){
  temp <- southern%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
  ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(1), Index=index.names[i])
  temp2<-rbind(ov3)
  overlap.southern <-rbind(temp2,overlap.southern)
}
overlap.southern


COP<-southern%>%bind_rows(northern)%>%mutate(yfirst=1989, ylast=2023,period=as.numeric(period))


####### Combining data for full plots #####
#full.dfa.spring<-RREAS%>%bind_rows(SEA,CALCOFI)

fulldat<-RREAS%>%bind_rows(CALCOFI, COP)%>%
  mutate(survey = fct_relevel(survey, "Northern Copepod (NCC)", 
                              "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"),
         period2 = ifelse(period==1,"1967 - 1988", ifelse(period==2, "1989 - 2012", "2013 - 2022")))

full.datslope<- ggplot(fulldat, aes(x = beta, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope",
       y = "Posterior density")

pdf(file = "Output/Figures/biological.slope.winter.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
full.datslope
dev.off()

full.datint<- ggplot(fulldat, aes(x = alpha, fill = as.factor(period2), group=as.factor(period2))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1], col[2], col[3]), name="Period") +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept",
       y = "Posterior density")

pdf(file = "Output/Figures/biological.winter.intercept.pdf",   # The directory you want to save the file in
    width = 8.5, # The width of the plot in inches
    height = 6)
full.datint
dev.off()


#### Violin plot data file #####

ratio.beta<-rbind(CALCOFI%>%filter(period=='3')%>%add_column((CALCOFI%>%filter(period=='3')%>%group_by(Index))$beta-(CALCOFI%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  northern%>%filter(period=='3')%>%add_column((northern%>%filter(period=='3')%>%group_by(Index))$beta-(northern%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  southern%>%filter(period=='3')%>%add_column((southern%>%filter(period=='3')%>%group_by(Index))$beta-(southern%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                  RREAS%>%filter(period=='3')%>%add_column((RREAS%>%filter(period=='3')%>%group_by(Index))$beta-(RREAS%>%filter(period=='2')%>%group_by(Index))$beta)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`))

colnames(ratio.beta)<- c("alpha", "beta", "survey","Index", "ratio")

ratio.alpha<-rbind(CALCOFI%>%filter(period=='3')%>%add_column((CALCOFI%>%filter(period=='3')%>%group_by(Index))$alpha-(CALCOFI%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                   northern%>%filter(period=='3')%>%add_column((northern%>%filter(period=='3')%>%group_by(Index))$alpha-(northern%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                   southern%>%filter(period=='3')%>%add_column((southern%>%filter(period=='3')%>%group_by(Index))$alpha-(southern%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`),
                   RREAS%>%filter(period=='3')%>%add_column((RREAS%>%filter(period=='3')%>%group_by(Index))$alpha-(RREAS%>%filter(period=='2')%>%group_by(Index))$alpha)%>%dplyr::select(alpha,beta,survey,Index,`... - ...`))

colnames(ratio.alpha)<- c("alpha", "beta", "survey", "Index","ratio")

# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }
ratio.alpha<-ratio.alpha%>%
  mutate(parameter="Intercept")

ratio.beta<-ratio.beta%>%
  mutate(parameter="Slope")

ratio<-ratio.beta%>%add_row(ratio.alpha)%>%
  mutate(survey2 = ifelse(survey=="CALCOFI (SCC)", "CalCOFI",
                          ifelse(survey=="RREAS (CCC)", "RREAS",
                                 ifelse(survey=="Southern Copepod (NCC)","s. Copepod","n. Copepod"))))%>%
  mutate(survey2 = factor(survey2, levels=c(c("CalCOFI","RREAS", "s. Copepod","n. Copepod"))))

col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

violin <-ggplot(ratio, aes(x=survey2, y=ratio, fill=survey2)) +
  theme_bw() +
  scale_fill_manual(values=col4, name="Survey")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~parameter,scales = "free_x") +
  ylab("Difference Era 3 - Era 2") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) 
pdf(file = "Output/Figures/violin.winter.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
violin
dev.off()
