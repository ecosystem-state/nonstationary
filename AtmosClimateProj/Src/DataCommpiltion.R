library(dplyr)
library(nord)
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
library(here)
library(PNWColors)
library(maps)       #basic mapping functions and some data
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(rgdal)
library(colorspace)
library(PBSmapping) #powerful mapping functions developed by Pacific Biological Station
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#### Importing Upwelling ####

here::i_am("Output/BayesianLinearModels.Rmd")


bakundat <-read.csv(here('data/physical/Bakun/erdUI246hr_d68d_e898_8529.csv'))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI276hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI306hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI336hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI366hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI396hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI426hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI456hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI486hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI516hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI546hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI576hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI606hr_d68d_e898_8529.csv')))%>%
  bind_rows(read.csv(here('data/physical/Bakun/erdUI616hr_d68d_e898_8529.csv')))


bakun <- bakundat%>%
  add_column('Year'=as.numeric(format(as.Date(bakundat$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(bakundat$time),"%m")))%>%
  add_column("Day"=as.numeric(format(as.Date(bakundat$time),"%d")))

bakun_region <- bakun%>%
  filter(station_id=='36N'|station_id=='39N')%>%
  mutate(region="Central CC")%>%
  bind_rows(bakun%>%
              filter(station_id=='33N'|station_id=='30N'|station_id=='27N'|station_id=='24N')%>%
              mutate(region="Southern CC"))%>%
  bind_rows(bakun%>%
              filter(station_id=='42N'|station_id=='45N'|station_id=='48N')%>%
              mutate(region="Northern CC"))%>%
  bind_rows(bakun%>%
              filter(station_id=='51N'|station_id=='54N'|station_id=='57N'|station_id=='60N')%>%
              mutate(region="GoA"))

bakun_season <- bakun_region%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(bakun_region%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(bakun_region%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))


bakun_summ_season2 <- bakun_season%>%  
  group_by(season, region) %>%
  summarise(seasonal_mean = mean(na.omit(upwelling_index)), seasonal_sd =sd(na.omit(upwelling_index)))

bakun_summ_annual2 <- bakun_season%>%  
  group_by(region) %>%
  summarise(yearly_mean = mean(na.omit(upwelling_index)), yearly_sd =sd(na.omit(upwelling_index)))

bakun_summ_annual <- bakun_season%>%
  group_by(region, Year_lag) %>%
  summarise(annual_mean = mean(na.omit(upwelling_index)))%>%
  ungroup()%>%
  right_join(bakun_summ_annual2, by = c( 'region'))%>%
  mutate(stand_bakun_annual = (annual_mean-yearly_mean)/yearly_sd)%>%
  select(Year_lag, region,stand_bakun_annual)

bakun_summ_seasonal <- bakun_season%>%
  group_by(region, Year_lag,  season) %>%
  summarise(season_mean = mean(na.omit(upwelling_index)))%>%
  ungroup()%>%
  right_join(bakun_summ_season2, by = c('season', 'region'))%>%
  mutate(stand_bakun_seasonally = (season_mean-seasonal_mean)/seasonal_sd)%>%
  select(season,Year_lag, region,stand_bakun_seasonally)

bakun_summ<- bakun_season%>%
  group_by(Month, region, Year_lag) %>%
  summarise(monthly_mean = mean(na.omit(upwelling_index)))%>%
  mutate(stand_bakun_annually = (monthly_mean-mean(na.omit(monthly_mean)))/sd(na.omit(monthly_mean)))%>%
  mutate(season=if_else(Month == 11|Month ==12|Month ==1|Month ==2|Month ==3, "Winter",
                        if_else(Month ==4|Month ==5|Month ==6, "Spring",
                                if_else(Month ==7|Month ==8, "Summer", "Autumn"))))%>%
  right_join(bakun_summ_season2, by = c('region', 'season'))%>%
  mutate(stand_bakun_monthly = (monthly_mean-seasonal_mean)/seasonal_sd)%>%
  left_join(bakun_summ_seasonal)%>%
  left_join(bakun_summ_annual)


bakun_time <-bakun_summ%>%
  filter(Year_lag>1963 & Year_lag<1989)%>%
  mutate(period='1')%>%
  bind_rows(bakun_summ%>%
              filter(Year_lag>1989 & Year_lag<2014)%>%
              mutate(period='2'))%>%
  bind_rows(bakun_summ%>%
              filter(Year_lag>2013)%>%
              mutate(period='3'))%>%
  mutate(era.region = if_else(region == "GoA"&period == 1, 1, 
                        if_else(region == "GoA"&period == 2, 2, 
                           if_else(region == "GoA"&period == 3, 3,
                              if_else(region == "Northern CC"&period == 1, 4,
                                  if_else(region == "Northern CC"&period == 2, 5,
                                      if_else(region == "Northern CC"&period == 3, 6,
                                           if_else(region == "Central CC"&period == 1, 7,
                                                if_else(region == "Central CC"&period == 2, 8,
                                                     if_else(region == "Central CC"&period == 3, 9, 
                                                          if_else(region == "Southern CC"&period == 1, 10,
                                                              if_else(region == "Southern CC"&period == 2, 11,
                                                                   if_else(region == "Southern CC"&period == 3, 12,
                                                                        13)))))))))))))
##### PDO ##### 
#### Climate Indices ####

PDO <- read.csv(here('data/physical/PDO.csv'))%>%
  pivot_longer(!Year, names_to = "Month", values_to = "PDO")%>%
  mutate(Month = as.integer(factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

PDO_seasonal <-  PDO%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(PDO%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(PDO%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_PDO = mean(PDO))


PDO_annual <-  PDO%>%
  group_by(Year_lag)%>%
  summarise(annual_PDO = mean(PDO))

##### NPGO #####
NPGO <- read.csv(here('data/physical/NPGO.csv'))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

NPGO_seasonal <-  NPGO%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(NPGO%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(NPGO%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_NPGO = mean(NPGO))


NPGO_annual <-  NPGO%>%
  group_by(Year_lag)%>%
  summarise(annual_NPGO = mean(NPGO))

##### ENSO #####

ONI <- read.csv(here('data/physical/ONI.csv'))%>%
  pivot_longer(!Year, names_to = "Month", values_to = "ONI")%>%
  mutate(Month = as.integer(factor(Month, levels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))))%>%
  filter(Year>1941)%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))

ONI_seasonal <-  ONI%>%
  filter(Month==12|Month==11|Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(ONI%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(ONI%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  mutate(Year_lag = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  group_by(Year_lag, season)%>%
  summarise(seasonal_ONI = mean(ONI))


ONI_annual <-  ONI%>%
  group_by(Year_lag)%>%
  summarise(annual_ONI = mean(ONI))

### Compile  with Upwelling Indices to Dataframe ####

climate_dat <- bakun_time%>%
  merge(PDO, by=c('Month', 'Year_lag'))%>%
  merge(PDO_annual, by=c('Year_lag'))%>%
  merge(PDO_seasonal, by=c('season', 'Year_lag'))%>%
  merge(ONI, by=c('Month', 'Year_lag'))%>%
  merge(ONI_annual, by=c('Year_lag'))%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO, by=c('Month', 'Year_lag'))%>%
  left_join(NPGO_annual, by=c('Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))


climate_dat <- climate_dat%>%
  select(Year_lag, region, season, seasonal_mean, stand_bakun_annual,period, 
         era.region,annual_PDO, seasonal_PDO, seasonal_NPGO, seasonal_ONI)%>%
  distinct()%>%
  filter(Year_lag<2023)

ggplot(climate_dat, aes(x=seasonal_PDO, y=seasonal_ONI))+
  geom_point()

saveRDS(climate_dat, file = here('data/physical/climate_dat_upwelling.rds'))

### Compile  Indices to Dataframe ####

climate_dat <-PDO_seasonal%>%
  merge(ONI_seasonal, by=c('season', 'Year_lag'))%>%
  left_join(NPGO_seasonal, by=c('season', 'Year_lag'))

climate_dat <-climate_dat%>%  
filter(Year_lag<1989)%>%
  mutate(period='1')%>%
  bind_rows(climate_dat%>%
              filter(Year_lag>1989 & Year_lag<2012)%>%
              mutate(period='2'))%>%
  bind_rows(climate_dat%>%
              filter(Year_lag>2012)%>%
              mutate(period='3'))

climate_dat <-climate_dat%>%  
  filter(Year_lag<1979)%>%
  mutate(period2='1')%>%
  bind_rows(climate_dat%>%
              filter(Year_lag>=1979)%>%
              mutate(period2='2'))%>%
  select(Year_lag, season, period, period2,seasonal_PDO, seasonal_NPGO, seasonal_ONI)%>%
  distinct()%>%
  filter(Year_lag<2023)

ggplot(climate_dat, aes(x=seasonal_PDO, y=seasonal_ONI))+
  geom_point()

saveRDS(climate_dat, file = here('data/physical/climate_dat.rds'))
