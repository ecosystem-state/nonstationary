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
bayeslinmod <- function(dat, ind,xyz, period) {
  N <- length(period)
  NP <- as.numeric(length(unique(period)))
  P<- as.numeric(as.factor(period))
  K <- 1
  y <-xyz
  x<-dat[ind]
  data <- list(N = N, #total number of opservations
               NP=NP, #total number of time periods
               P = P, #time period pointer vector
               K = K,#number of covariates, starting with one but can add for final model structure
               x=x, 
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
  alphaP<-posterior[1:NP,]%>%
    add_column(
      period=rep(seq(1,NP),1))
  
  betaP<-posterior[NP+1:NP*2,]%>%
    add_column(
      period=rep(seq(1,NP),1))
  n<- 10000
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled <-rbind(scaled,cbind(tempalpha, tempbeta))
  }
  scaled.anom<<- scaled%>%
    bind_cols(period=rep(rep(seq(1,NP), each = n),1), Index=rep(colnames(dat[ind]), n*NP),
              yfirst=rep(min(dat$Year_lag), n*NP),ylast=rep(max(dat$Year_lag), n*NP))%>%
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
#### CALCOFI ####
####Spring ####

data <- dat%>%filter(season=="Spring", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

CALCOFI <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  CALCOFI <-rbind(CALCOFI,scaled.anom)
}

CALCOFI<-CALCOFI%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI, CALCOFI$beta)
postplot(CALCOFI, CALCOFI$alpha)

CALCOFI_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_full <-rbind(CALCOFI_full,scaled.anom)
}

CALCOFI_full<-CALCOFI_full%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI<- bind_rows(CALCOFI,CALCOFI_full)

index.names <- unique(CALCOFI$Index)
period.names <- unique(CALCOFI$period)
overlap.CALCOFI <- NA
  for(i in 1:4){
    temp <- CALCOFI%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFI <-rbind(temp2,overlap.CALCOFI)
  }
overlap.CALCOFI 

#running models with 1-year offset

CALCOFIoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  CALCOFIoffset <-rbind(CALCOFIoffset,scaled.anom)
}
CALCOFIoffset<-CALCOFIoffset%>%mutate(survey="CALCOFI (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(CALCOFIoffset, CALCOFIoffset$beta)
postplot(CALCOFIoffset, CALCOFIoffset$alpha)

CALCOFIoffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFIoffset_full <-rbind(CALCOFIoffset_full,scaled.anom)
}

CALCOFIoffset_full<-CALCOFIoffset_full%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFIoffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFIoffset<- bind_rows(CALCOFIoffset,CALCOFIoffset_full)

index.names <- unique(CALCOFIoffset$Index)
period.names <- unique(CALCOFIoffset$period)
overlap.CALCOFIoffset <- NA
  for(i in 1:4){
    temp <- CALCOFIoffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFIoffset <-rbind(temp2,overlap.CALCOFIoffset)
  }
overlap.CALCOFIoffset 

#### Winter ####

data <- dat%>%filter(season=="Winter", trend=="CALCOFI")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

CALCOFI_W <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  CALCOFI_W <-rbind(CALCOFI_W,scaled.anom)
}

CALCOFI_W<-CALCOFI_W%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(CALCOFI_W, CALCOFI_W$beta)
postplot(CALCOFI_W, CALCOFI_W$alpha)

CALCOFI_W_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_W_full <-rbind(CALCOFI_W_full,scaled.anom)
}

CALCOFI_W_full<-CALCOFI_W_full%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_W_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI_W<- bind_rows(CALCOFI_W,CALCOFI_W_full)

index.names <- unique(CALCOFI_W$Index)
period.names <- unique(CALCOFI_W$period)
overlap.CALCOFI_W <- NA
  for(i in 1:4){
    temp <- CALCOFI_W%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFI_W <-rbind(temp2,overlap.CALCOFI_W)
  }
overlap.CALCOFI_W 

#running models with 1-year offset

CALCOFI_Woffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  CALCOFI_Woffset <-rbind(CALCOFI_Woffset,scaled.anom)
}
CALCOFI_Woffset<-CALCOFI_Woffset%>%mutate(survey="CALCOFI (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(CALCOFI_Woffset, CALCOFI_Woffset$beta)
postplot(CALCOFI_Woffset, CALCOFI_Woffset$alpha)

CALCOFI_Woffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  CALCOFI_Woffset_full <-rbind(CALCOFI_Woffset_full,scaled.anom)
}

CALCOFI_Woffset_full<-CALCOFI_Woffset_full%>%mutate(survey="CALCOFI (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(CALCOFI_Woffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


CALCOFI_Woffset<- bind_rows(CALCOFI_Woffset,CALCOFI_Woffset_full)

index.names <- unique(CALCOFI_Woffset$Index)
period.names <- unique(CALCOFI_Woffset$period)
overlap.CALCOFI_Woffset <- NA
  for(i in 1:4){
    temp <- CALCOFI_Woffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.CALCOFI_Woffset <-rbind(temp2,overlap.CALCOFI_Woffset)
  }
overlap.CALCOFI_Woffset 



#### RREAS ####
#####Spring #####

data <- dat%>%filter(season=="Spring", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

RREAS <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  RREAS <-rbind(RREAS,scaled.anom)
}

RREAS<-RREAS%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(RREAS, RREAS$beta)
postplot(RREAS, RREAS$alpha)

RREAS_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_full <-rbind(RREAS_full,scaled.anom)
}

RREAS_full<-RREAS_full%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS<- bind_rows(RREAS,RREAS_full)

index.names <- unique(RREAS$Index)
period.names <- unique(RREAS$period)
overlap.RREAS <- NA
  for(i in 1:4){
    temp <- RREAS%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS <-rbind(temp2,overlap.RREAS)
  }
overlap.RREAS 

#running models with 1-year offset

RREASoffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  RREASoffset <-rbind(RREASoffset,scaled.anom)
}
RREASoffset<-RREASoffset%>%mutate(survey="RREAS (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(RREASoffset, RREASoffset$beta)
postplot(RREASoffset, RREASoffset$alpha)

RREASoffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREASoffset_full <-rbind(RREASoffset_full,scaled.anom)
}

RREASoffset_full<-RREASoffset_full%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREASoffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREASoffset<- bind_rows(RREASoffset,RREASoffset_full)

index.names <- unique(RREASoffset$Index)
period.names <- unique(RREASoffset$period)
overlap.RREASoffset <- NA
  for(i in 1:4){
    temp <- RREASoffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREASoffset <-rbind(temp2,overlap.RREASoffset)
  }
overlap.RREASoffset 

#### Winter ####

data <- dat%>%filter(season=="Winter", trend=="RREAS")%>%
  distinct()%>%
  filter(Year_lag<2023)%>%
  add_column(full=1)

RREAS_W <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$period)
  RREAS_W <-rbind(RREAS_W,scaled.anom)
}

RREAS_W<-RREAS_W%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

postplot(RREAS_W, RREAS_W$beta)
postplot(RREAS_W, RREAS_W$alpha)

RREAS_W_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_W_full <-rbind(RREAS_W_full,scaled.anom)
}

RREAS_W_full<-RREAS_W_full%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS_W_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS_W<- bind_rows(RREAS_W,RREAS_W_full)

index.names <- unique(RREAS_W$Index)
period.names <- unique(RREAS_W$period)
overlap.RREAS_W <- NA
  for(i in 1:4){
    temp <- RREAS_W%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS_W <-rbind(temp2,overlap.RREAS_W)
  }
overlap.RREAS_W 

#running models with 1-year offset

RREAS_Woffset <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimateoffset1, data$period)
  RREAS_Woffset <-rbind(RREAS_Woffset,scaled.anom)
}
RREAS_Woffset<-RREAS_Woffset%>%mutate(survey="RREAS (SCC)",
                          Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                            ifelse(Index=="seasonal_NPH", "NPH","ONI"))))
  
postplot(RREAS_Woffset, RREAS_Woffset$beta)
postplot(RREAS_Woffset, RREAS_Woffset$alpha)

RREAS_Woffset_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$estimate, data$full)
  RREAS_Woffset_full <-rbind(RREAS_Woffset_full,scaled.anom)
}

RREAS_Woffset_full<-RREAS_Woffset_full%>%mutate(survey="RREAS (SCC)",
Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                  ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  mutate(period=4)

ggplot(RREAS_Woffset_full, aes(x = beta),fill = as.factor(period), group=as.factor(period)) +
    theme_bw() +
    facet_wrap(.~Index, ncol = 4, scales='free') +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(col[1])) +
    #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
    geom_vline(xintercept = 0, lty = 2) +
    labs(x = "Slope",
         y = "Posterior density")


RREAS_Woffset<- bind_rows(RREAS_Woffset,RREAS_Woffset_full)

index.names <- unique(RREAS_Woffset$Index)
period.names <- unique(RREAS_Woffset$period)
overlap.RREAS_Woffset <- NA
  for(i in 1:4){
    temp <- RREAS_Woffset%>%filter(Index==index.names[i])%>%dplyr::select(beta, period)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(4), period2=c(3), Index=index.names[i])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==1)%>%dplyr::select(beta)), period1=c(4), period2=c(2), Index=index.names[i])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==4)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(4), period2=c(1), Index=index.names[i])

    temp2<-rbind(ov1,ov2,ov3, ov4, ov5, ov6)
    overlap.RREAS_Woffset <-rbind(temp2,overlap.RREAS_Woffset)
  }
overlap.RREAS_Woffset 


