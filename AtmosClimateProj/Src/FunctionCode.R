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
climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

dat.long<- dat%>%filter(season=="Spring")%>%
  select(Year_lag, season, estimate,trend,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

dat.lm<-climate_dat_cop%>%filter(season=="Spring")%>%
  select(Year_lag, season, seasonal_copepod_northern,seasonal_copepod_southern,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag,  season,NPH,NPGO,PDO,ONI), 
               values_to = "estimate", names_to = "trend")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
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
season<- "Spring"
data <- dat%>%filter(season==season, trend=="CALCOFI")%>%
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

data <- dat%>%filter(season==season, trend=="RREAS")%>%
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
  n<- 1000
  scaled.anomaly <- NULL
  for(i in 1:NP){
    tempalpha <- rnorm(n, alphaP$mean[i],alphaP$sd[i])
    tempbeta <- rnorm(n, betaP$mean[i],betaP$sd[i])
    scaled.anomaly <-rbind(scaled.anomaly,cbind(tempalpha, tempbeta))
  }
  
  scaled.anom<<- scaled.anomaly%>%
    bind_cols(period=rep(rep(c('2', '3'), each = 1000),1), Index = rep(colnames(dat[ind]), 2000))%>%
    rename(alpha=tempalpha, beta = tempbeta)
}


climpdo <- climate_dat_cop%>%filter(season=="Spring")%>%
  select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
  distinct()%>%
  filter(Year_lag<2023)

data <- climate_dat_cop%>%filter(season=="Summer")%>%
  select(Year_lag, region, seasonal_copepod_northern, annual_copepod_northern, 
         seasonal_copepod_southern)%>%
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
northern<-northern%>%mutate(survey="Northern Copepod (NCC)",
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                                              ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

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
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
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
up_post_spring<-up_post%>%
  mutate(era=as.numeric(as.factor(period)),
         period2 = ifelse(era==1,"1967 - 1988", ifelse(era==2, "1989 - 2012", "2013 - 2023")),
         Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                                           ifelse(Index=="seasonal_NPH", "NPH","ONI"))))

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


dat.lm<- data%>%select(Year_lag, period,region, season, era.region,stand_bakun_seasonally,seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
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
