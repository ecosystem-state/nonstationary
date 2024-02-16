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

dat<- readRDS("data/all_juvenile_indices.rds")
dat<- read.csv("data/salmon/MatRate_Survival_up to CY 2021_new.csv")

Indicator_ID<-data.frame(stock_code=unique(dat$stock_code))
ID1<-Indicator_ID%>%mutate(stock_name=ifelse(stock_code=='ATN', "Atnarko", 
                                 ifelse(stock_code=='BQR', "Big Qualicum River Fall",
                                 ifelse(stock_code=='CHI', "Chilliwack River Fall",
                                 ifelse(stock_code=='COW', "Cowichan River Fall",
                                 ifelse(stock_code=='CWF', "Cowlitz Fall Tule",
                                 ifelse(stock_code=='ELK', "Elk River",
                                 ifelse(stock_code=='ELW', "Elwha River",
                                 ifelse(stock_code=='GAD', "George Adams Fall Fingerling",
                                 ifelse(stock_code=='GRN', "Green River Fingerling",
                                 ifelse(stock_code=='HAN', "Hanford Wild Brights",
                                 ifelse(stock_code=='HAR', "Harrison River",
                                 ifelse(stock_code=='HOK', "Hoko Fall Fingerling",
                                 ifelse(stock_code=='KLM', "Kitsumkalum River Summer",
                                 ifelse(stock_code=='KLY', "Kitsumkalum Yearling",
                                 ifelse(stock_code=='LRH', "Columbia Lower River Hatchery",
                                 ifelse(stock_code=='LRW', "Lewis River Wild",
                                 ifelse(stock_code=='LYF', "Lyons Ferry Fingerling",
                                 ifelse(stock_code=='LYY', "Lyons Ferry Yearling", 
                                 ifelse(stock_code=='MSH', "Middle Shuswap",
                                 ifelse(stock_code=='NIC', "Nicola",
                                 ifelse(stock_code=='NIS', "Nisqually Fall Fingerling",
                                 ifelse(stock_code=='NSA', "Northern Southeast Alaska",
                                 ifelse(stock_code=='NSF', "Nooksack Spring Fingerling",
                                 ifelse(stock_code=='PHI', "Phillips River Fall",
                                 ifelse(stock_code=='PPS', "Puntledge River Summer",
                                 ifelse(stock_code=='QUE', "Queets Fall Fingerling",
                                 ifelse(stock_code=='QUI', "Quinsam River Fall ",
                                 ifelse(stock_code=='RBT', "Robertson Creek Fall",
                                 ifelse(stock_code=='SAM', "Samish Fall Fingerling",
                                 ifelse(stock_code=='SHU', "Lower Shuswap River Summer",
                                 ifelse(stock_code=='SKF', "Skagit Spring Fingerling",
                                 ifelse(stock_code=='SKY', "Skykomish Fall Fingerling",
                                 ifelse(stock_code=='SMK', "Similkameen Summer Yearling",
                                 ifelse(stock_code=='SOO', "Tsoo-Yess Fall Fingerling",
                                 ifelse(stock_code=='SPR', "Spring Creek Tule",
                                 ifelse(stock_code=='SPS', "South Puget Sound Fall Fingerling",
                                 ifelse(stock_code=='SRH', "Salmon River",
                                 ifelse(stock_code=='SSA', "Southern SEAK Spring",
                                 ifelse(stock_code=='SSF', "Skagit Summer Fingerling",
                                 ifelse(stock_code=='STL', "Stillaguamish Fall Fingerling² ",
                                 ifelse(stock_code=='SUM', "Columbia River Summers",
                                 ifelse(stock_code=='URB', "Columbia Upriver Bright",
                                 ifelse(stock_code=='WSH', "Willamette Spring",1))))))))))))))))))))))))))))))))))))))))))))


ID2<-data.frame(Indicator_ID%>%mutate(area=ifelse(stock_code=='NSA'|stock_code=='SSA'|stock_code=='CHK'|stock_code=='UNU', "Southeast Alaska", 
                                 ifelse(stock_code=='TST', "Transboundary Rivers",
                                 ifelse(stock_code=='ATN'|stock_code=='KLM'|stock_code=='KLY', "North/Central BC",
                                 ifelse(stock_code=='RBT', "WCVI",
                                 ifelse(stock_code=='QUI'|stock_code=='PHI'|stock_code=='PPS'|stock_code=='BQR'|stock_code=='COW', "Strait of Georgia",
                                 ifelse(stock_code=='HAR'|stock_code=='CHI'|stock_code=='CKO'|stock_code=='NIC'|stock_code=='SHU'|stock_code=='MSH', "Fraser River",
                                 ifelse(stock_code=='NSF'|stock_code=='SAM'|stock_code=='SSF'|stock_code=='SKF', "North Puget Sound",
                                 ifelse(stock_code=='STL'|stock_code=='SKY', "Central Puget Sound",
                                 ifelse(stock_code=='NIS'|stock_code=='SPS'|stock_code=='GRN', "South Puget Sound",
                                 ifelse(stock_code=='GAD', "Hood Canal",
                                 ifelse(stock_code=='ELW', "Juan de Fuca",
                                 ifelse(stock_code=='HOK'|stock_code=='QUE'|stock_code=='SOO', "North Washington Coast",
                                 ifelse(stock_code=='LRH'|stock_code=='CWF'|stock_code=='LRW'|stock_code=='WSH'|stock_code=='SPR', "Lower Columbia River",
                                 ifelse(stock_code=='HAN'|stock_code=='SMK'|stock_code=='SUM'|stock_code=='URB', "Upper Columbia River",
                                 ifelse(stock_code=='LYF'|stock_code=='LYY', "Snake River",
                                 ifelse(stock_code=='SRH', "North Oregon Coast",
                                 ifelse(stock_code=='ELK', "Mid Oregon Coast",1)))))))))))))))))))

ID3<-ID2%>%mutate(ecoregion=ifelse(area=="North/Central BC", 2, 
                                 ifelse(area=="Strait of Georgia", 2,
                                 ifelse(area=="Fraser River", 2,
                                 ifelse(area=="Lower Columbia River", 2,
                                 ifelse(area=="Mid Oregon Coast", 3,
                                 ifelse(area=="Juan de Fuca", 2,
                                 ifelse(area=="Hood Canal", 1,
                                 ifelse(area=="South Puget Sound", 1,
                                 ifelse(area=="Upper Columbia River", 2,
                                 ifelse(area=="North Washington Coast", 2,
                                 ifelse(area=="Snake River", 2,
                                 ifelse(area=="Southeast Alaska", 0,
                                 ifelse(area=="North Puget Sound",1,
                                 ifelse(area=="WCVI", 0,
                                 ifelse(area=="Central Puget Sound", 1,
                                 ifelse(area=="North Oregon Coast", 3,0)))))))))))))))))

Marine_Survival<-dat%>%filter(age==2)%>%mutate(calendar_year=brood_year+1)%>%
  left_join(ID1)%>%left_join(ID2)%>%left_join(ID3)%>%mutate(location=ifelse(ecoregion==1|ecoregion==2,"47N",
                                                            ifelse(ecoregion==3,"45N",NA)))%>%
  mutate(location2="45N")

  


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

Beuti<-read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_45N.csv')%>%add_column(location2='45N')%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_38N.csv')%>%add_column(location2='38N'))%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_BEUTI_33N.csv')%>%add_column(location2='33N'))

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
  group_by(Year, season, location2)%>%
  summarise(seasonal_beuti = mean(beuti))%>%
  ungroup()%>%
  group_by(season, location2)%>%
  mutate(seasonal_beuti=scale(seasonal_beuti))%>%
  rename(calendar_year=Year)%>%
  ungroup()


Beuti_seasonal%>%summarise(mean=mean(seasonal_beuti), sd=sd(seasonal_beuti))

CUTI<-read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_45N.csv')%>%add_column(location2='45N')%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_38N.csv')%>%add_column(location2='39N'))%>%
  bind_rows(read.csv('data/Environment/CUTI_BEUTI/cciea_OC_CUTI_33N.csv')%>%add_column(location2='33N'))

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
  group_by(Year, season, location2)%>%
  summarise(seasonal_cuti = mean(cuti))%>%
  ungroup()%>%
  group_by(season, location2)%>%
  mutate(seasonal_cuti=scale(seasonal_cuti))%>%
  rename(calendar_year=Year)%>%
  ungroup()

STI<-read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_39N.csv')%>%add_column(location2='39N')%>%
  bind_rows(read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_42N.csv')%>%add_column(location2='42N'))%>%
  bind_rows(read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_45N.csv')%>%add_column(location2='45N'))%>%
  bind_rows(read.csv('data/Environment/Upwelling_Phenology/cciea_OC_STI_48N.csv')%>%add_column(location2='48N'))

STI<-STI%>%
  add_column('Year'=as.numeric(format(as.Date(STI$time),"%Y")))%>%
  add_column('Month'=as.numeric(format(as.Date(STI$time),"%m")))%>%
  group_by(Year, location2)%>%
  summarise(seasonal_STI = mean(sti))%>%
  ungroup()%>%
  group_by(location2)%>%
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
  dplyr::left_join(CUTI_winter, by=c('location2', 'calendar_year'))%>%
  dplyr::left_join(CUTI_spring, by=c('location2', 'calendar_year'))%>%
 dplyr::left_join(STI, by=c('location2', 'calendar_year'))%>%
  dplyr::left_join(Beuti_winter, by=c('location2', 'calendar_year'))%>%
  dplyr::left_join(Beuti_spring, by=c('location2', 'calendar_year'))%>%
  dplyr::left_join(Beuti_eco_winter, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Beuti_eco_spring, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Bifurcation, by=c('calendar_year'))
write.csv(Survival_combined, file = "data/Survival_combined.csv")

##### Exploratory Plots #####
data<-read.csv("data/Survival_combined.csv")

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
