#### Upwelling Model Runs ####

##### spring #####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
season <- "Spring"
data<- climate_dat%>%filter(season=="Spring")%>%
  #  filter(Year_lag!=2022&Year_lag!=2015)%>%
  distinct()

upwelling <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$era.region)
  upwelling <-rbind(upwelling,scaled.anom)
}

upwelling<-upwelling%>%mutate(survey="Upwelling",era.region=period,
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  dplyr::select(!period)%>%
  left_join(data%>%dplyr::select(region,era.region,period)%>%distinct())%>%
  mutate(Season='Spring', lag=0)


unique(upwelling$period)



upwelling_full <-NULL
columns<-c(which(colnames(data) == "seasonal_PDO"), which(colnames(data) == "seasonal_NPGO"),
           which(colnames(data) == "seasonal_NPH"),which(colnames(data) == "seasonal_ONI"))
for(i in 1:length(columns)){
  bayeslinmod(data, columns[i],  data$stand_bakun_seasonally, data$region)
  upwelling_full <-rbind(upwelling_full,scaled.anom)
}

reg<-unique(data%>%dplyr::select(region))%>%
  rename(region2=region)%>%
  mutate(region=as.numeric(as.factor(region2)))

upwelling_full<-upwelling_full%>%mutate(survey="Upwelling",region=period,
                            Index=ifelse(Index=="seasonal_PDO", "PDO", ifelse(Index=="seasonal_NPGO", "NPGO",
                                ifelse(Index=="seasonal_NPH", "NPH","ONI"))))%>%
  left_join(reg)%>%mutate(period=4)%>%dplyr::select(!region)%>%rename(region=region2)%>%
  mutate(Season='Spring', lag=0)


upwelling <- upwelling%>%dplyr::select(!era.region)%>%mutate(period=as.numeric(as.factor(period)))%>%
  bind_rows(upwelling_full%>%dplyr::select(alpha,beta,Index,yfirst,ylast,survey,region,period,Season,lag))

index.names <- unique(upwelling$Index)
region.names<-unique(upwelling$region)
period.names <- unique(upwelling$period)
overlap.up <- NA
for(j in 1:4){
  temp1 <- upwelling%>%filter(region==region.names[j])
  for(i in 1:4){
    temp <- temp1%>%filter(Index==index.names[i])%>%dplyr::select(beta, period,region)
    ov1 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(1), period2=c(2), Index=index.names[i], region=region.names[j])
    ov2 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==3)%>%dplyr::select(beta)), period1=c(1), period2=c(3), Index=index.names[i], region=region.names[j])
    ov3 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==2)%>%dplyr::select(beta)), period1=c(3), period2=c(2), Index=index.names[i], region=region.names[j])
    ov4 <- data.frame(ov=overlap(temp%>%filter(period==1)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(1), period2=c(4), Index=index.names[i], region=region.names[j])
    ov5 <- data.frame(ov=overlap(temp%>%filter(period==2)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(2), period2=c(4), Index=index.names[i], region=region.names[j])
    ov6 <- data.frame(ov=overlap(temp%>%filter(period==3)%>%dplyr::select(beta),temp%>%filter(period==4)%>%dplyr::select(beta)), period1=c(3), period2=c(4), Index=index.names[i], region=region.names[j])

    temp2<-rbind(ov1,ov2,ov3,ov4,ov5,ov6)
    overlap.up <-rbind(temp2,overlap.up)
  }
  
}


