library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(mgcv)
library(MASS)
library(stringr)
library(gamm4)
library(tidyr)
library(ggthemes)
library(viridis)
library(cowplot)
library(kableExtra)
library(docxtools)
library(knitr)
library(tibble)
library(gratia)
library(latex2exp)
library(patchwork)

#dat<- readRDS("data/all_juvenile_indices.rds")
dat<- read.csv("data/salmon/MatRate_Survival_up to CY 2021_new.csv")%>%
  dplyr::select(stock_code, brood_year, age, Mat.Rate, Marine.Survival)
cwt<- read.csv("data/CWTinfo.csv")
Klamath<- read.csv("data/salmon/SurvAge2.csv")%>%
  rename(Marine.Survival=s.rate, brood_year=brood.yr, stock_code=component)%>%
  mutate(Mat.Rate=0)%>%
  dplyr::select(stock_code, brood_year, age, Mat.Rate, Marine.Survival)

Indicator_ID<-data.frame(stock_code=unique(dat$stock_code))
dat<-dat%>%
  add_row(Klamath)%>%
  left_join(cwt)
Marine_Survival<-dat%>%filter(age==2)%>%mutate(calendar_year=brood_year+1)%>%
  left_join(ID1)%>%left_join(ID2)%>%left_join(ID3)%>%left_join(ID4)%>%mutate(location=ifelse(ecoregion==1|ecoregion==2,"47N",
                                                            ifelse(ecoregion==3,"45N",
                                                                  ifelse(ecoregion==4,"39N", NA))))%>%
 mutate(location=ifelse(ecoregion==1|ecoregion==2,"45N",
                                                            ifelse(ecoregion==3,"45N",
                                                                  ifelse(ecoregion==4,"39N", NA))))



##### Environmental Data ####
Bifurcation<-read.csv('data/Environment/BifurcationIndex.csv')%>%
mutate(stand_BI=scale(BI))%>%
  rename(calendar_year=Year)

Beuti_eco<-read.csv('data/Environment/CUTI_BEUTI/BEUTI_monthly.csv')%>%
  dplyr::select(year, month,X47N, X45N, X42N,X39N)
Beuti_eco_seasonal<-Beuti_eco%>%
  filter(month==1|month==2|month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(Beuti_eco%>%
              filter(month==4|month==5|month==6)%>%
              mutate(season="Spring"))%>%
  pivot_longer(c(X47N,X45N,X42N,X39N),names_to = "location1",values_to="Beuti_eco" )%>%
  mutate(location = ifelse(location1 == 'X47N', "47N",
                    ifelse(location1 == 'X45N', "45N",
                    ifelse(location1 == 'X42N',"42N",
                    ifelse(location1 == 'X39N', "39N",
                         location1)))))%>%
  dplyr::select(-location1)%>%
  group_by(year, season, location)%>%
  summarise(seasonal_Beuti_eco = mean(Beuti_eco))%>%
  ungroup()%>%
  group_by(season, location)%>%
  mutate(seasonal_Beuti_eco=scale(seasonal_Beuti_eco))%>%
  rename(calendar_year=year)%>%
  ungroup()

Beuti_eco_seasonal%>%group_by(season, location)%>%
  summarise(mean=mean(seasonal_Beuti_eco), sd=sd(seasonal_Beuti_eco))

Beuti<-read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_45N.csv')%>%add_column(location='45N')%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_38N.csv')%>%add_column(location='39N'))%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_33N.csv')%>%add_column(location='33N'))

Beuti<-Beuti%>%
  add_column('Year'=as.numeric(format(as.Date(Beuti$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(Beuti$time),"%m")))
  
Beuti_seasonal<-Beuti%>%
 filter(Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(Beuti%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(Beuti%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  group_by(Year, season, location)%>%
  summarise(seasonal_beuti = mean(beuti))%>%
  ungroup()%>%
  group_by(season, location)%>%
  mutate(seasonal_beuti=scale(seasonal_beuti))%>%
  rename(calendar_year=Year)%>%
  ungroup()


Beuti_seasonal%>%summarise(mean=mean(seasonal_beuti), sd=sd(seasonal_beuti))

CUTI<-read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_45N.csv')%>%add_column(location='45N')%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_38N.csv')%>%add_column(location='39N'))%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_33N.csv')%>%add_column(location='33N'))

CUTI<-CUTI%>%
  add_column('Year'=as.numeric(format(as.Date(CUTI$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(CUTI$time),"%m")))

CUTI_seasonal<-CUTI%>%
 filter(Month==1|Month==2|Month==3)%>%
  mutate(season="Winter")%>%
  bind_rows(CUTI%>%
              filter(Month==4|Month==5|Month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(CUTI%>%
              filter(Month==7|Month==8)%>%
              mutate(season="Summer"))%>%
  group_by(Year, season, location)%>%
  summarise(seasonal_cuti = mean(cuti))%>%
  ungroup()%>%
  group_by(season, location)%>%
  mutate(seasonal_cuti=scale(seasonal_cuti))%>%
  rename(calendar_year=Year)%>%
  ungroup()

STI<-read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_39N.csv')%>%add_column(location='39N')%>%
  bind_rows(read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_42N.csv')%>%add_column(location='42N'))%>%
  bind_rows(read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_45N.csv')%>%add_column(location='45N'))

STI<-STI%>%
  add_column('Year'=as.numeric(format(as.Date(STI$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(STI$time),"%m")))%>%
  group_by(Year, location)%>%
  summarise(seasonal_STI = mean(sti))%>%
  ungroup()%>%
  group_by(location)%>%
  mutate(seasonal_STI=scale(seasonal_STI))%>%
  rename(calendar_year=Year)%>%
  ungroup()

SST.dat <-readRDS('data/environment/SST/ERSST_poly.rds')

SST_seasonal<-SST.dat%>%
 filter(month==1|month==2|month==3)%>%
  filter(year>=1975)%>%
  mutate(season="Winter")%>%
  bind_rows(SST.dat%>%
              filter(month==4|month==5|month==6)%>%
              mutate(season="Spring"))%>%
  bind_rows(SST.dat%>%
              filter(month==7|month==8)%>%
              mutate(season="Summer"))%>%
  group_by(year, season, ecoregion)%>%
  summarise(seasonal_sst = mean(sst))%>%
  ungroup()%>%
  group_by(season, ecoregion)%>%
  mutate(seasonal_sst=scale(seasonal_sst))%>%
  rename(calendar_year=year)%>%
  ungroup()

#### Combing Environmental with Chinook####
SST_winter<-SST_seasonal%>%filter(season=='Winter')%>%rename(sst_winter=seasonal_sst)%>%dplyr::select(-season)
SST_spring<-SST_seasonal%>%filter(season=='Spring')%>%rename(sst_spring=seasonal_sst)%>%dplyr::select(-season)

CUTI_winter<-CUTI_seasonal%>%filter(season=='Winter')%>%rename(CUTI_winter=seasonal_cuti)%>%dplyr::select(-season)
CUTI_spring<-CUTI_seasonal%>%filter(season=='Spring')%>%rename(CUTI_spring=seasonal_cuti)%>%dplyr::select(-season)

Beuti_winter<-Beuti_seasonal%>%filter(season=='Winter')%>%rename(Beuti_winter=seasonal_beuti)%>%dplyr::select(-season)
Beuti_spring<-Beuti_seasonal%>%filter(season=='Spring')%>%rename(Beuti_spring=seasonal_beuti)%>%dplyr::select(-season)

Beuti_eco_winter<-Beuti_eco_seasonal%>%filter(season=='Winter')%>%rename(Beuti_eco_winter=seasonal_Beuti_eco)%>%dplyr::select(-season)
Beuti_eco_spring<-Beuti_eco_seasonal%>%filter(season=='Spring')%>%rename(Beuti_eco_spring=seasonal_Beuti_eco)%>%dplyr::select(-season)


Survival_combined<-Marine_Survival%>%dplyr::left_join(SST_winter,by=c('ecoregion', 'calendar_year'))%>%
                                dplyr::left_join(SST_spring, by=c('ecoregion', 'calendar_year'))%>%
  dplyr::left_join(CUTI_winter, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(CUTI_spring, by=c('location', 'calendar_year'))%>%
 dplyr::left_join(STI, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Beuti_winter, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Beuti_spring, by=c('location', 'calendar_year'))%>%
 # dplyr::left_join(Beuti_eco_winter, by=c('location', 'calendar_year'))%>%
 # dplyr::left_join(Beuti_eco_spring, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Bifurcation, by=c('calendar_year'))
write.csv(Survival_combined, file = "data/Salmon/Survival_combined.csv")

##### Exploratory Plots #####
data<-read.csv("data/Survival_combined.csv")
unique(data$stock_code)
TS.survival<-ggplot(data =Survival_combined,
              aes(x =calendar_year, y = Marine.Survival, group = stock_name, col=stock_name))+
  #group=region)) +
  facet_wrap(.~area, ncol = 3, labeller = label_wrap_gen(25), scales="free_y") +
  geom_line(aes(group = stock_name,col=stock_name))+
  scale_y_continuous(name ="Age-2 Survival" )+
  scale_x_continuous(name = "Calendar Year")

pdf(file = "Results/Figures/Age2SurvivalTS.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 7)
TS.survival
dev.off()
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  his <- hist(x, plot = FALSE)
  breaks <- his$breaks
  nB <- length(breaks)
  y <- his$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = rgb(0, 1, 1, alpha = 0.5), ...)
  # lines(density(x), col = 2, lwd = 2) # Uncomment to add density lines
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  Cor <- abs(cor(x, y)) # Remove abs function if desired
  txt <- paste0(prefix, format(c(Cor, 0.123456789), digits = digits)[1])
  if(missing(cex.cor)) {
    cex.cor <- 0.4 / strwidth(txt)
  }
  text(0.5, 0.5, txt,
       cex = 1 + cex.cor * Cor) # Resize the text by level of correlation
}

pdf(file = "Results/Figures/CovariateCorrelation.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
pairs(~ seasonal_STI+
        Beuti_spring+
        stand_BI+
        CUTI_spring+
        sst_spring+
        Beuti_winter+
        CUTI_winter+
        sst_winter, 
      data = dd, upper.panel = panel.cor,         # Disabling the upper panel
      diag.panel = panel.hist,
      lower.panel = panel.smooth)

dev.off()
vif.covs<-car::vif(lm(Marine.Survival~seasonal_STI+Beuti_spring+stand_BI+sst_spring+Beuti_winter, data=dd))
write.table(vif.covs, file = "Results/Tables/VIF.txt", sep = ",", quote = FALSE, row.names = F)
