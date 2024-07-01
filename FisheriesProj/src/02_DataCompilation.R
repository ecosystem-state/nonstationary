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
library(here)

#### Reading in the salmon data ####
#dat<- readRDS("data/all_juvenile_indices.rds")
dat<- read.csv("data/salmon/MatRate_Survival_up to CY 2021_new.csv")%>% # data from CTC
  dplyr::select(stock_code, brood_year, age, Mat.Rate, Marine.Survival)
cwt<- read.csv("data/CWTinfo.csv")

Klamath<- read.csv("data/salmon/SurvAge2.csv")%>% #Klamath data from Will
  rename(Marine.Survival=s.rate, brood_year=brood.yr, stock_code=component)%>%
  mutate(Mat.Rate=NA)%>% #filler mat rate because we dont have the data for this stock
  dplyr::select(stock_code, brood_year, age, Mat.Rate, Marine.Survival)

Indicator_ID<-data.frame(stock_code=unique(dat$stock_code))
dat<-dat%>% #Combining Klamath and CTC into a single dataframe
  add_row(Klamath)%>%
  left_join(cwt)


#pulling out specifically age 2 survival and also assigning a location to align with 
#covariate data 'ecoregion' and also making calendar year = brood year + lag to align with the 
#covariate data
Marine_Survival<-dat%>%
  filter(age==2)%>%
  mutate(calendar_year=brood_year+1)%>%
#  left_join(ID1)%>%
#  left_join(ID2)%>%
#  left_join(ID3)%>%
#  left_join(ID4)%>%
  mutate(location=ifelse(ecoregion==1|ecoregion==2,"47N",
                  ifelse(ecoregion==3,"45N",
                  ifelse(ecoregion==4,"39N", NA))))%>%
 mutate(location=ifelse(ecoregion==1|ecoregion==2,"45N",
                 ifelse(ecoregion==3,"45N",
                  ifelse(ecoregion==4,"39N", NA))))%>%
  mutate(ecoregion=ifelse(ecoregion==1,2,ecoregion)) #pull this out for original ecoregion characterization


##### Accessing and summarising Environmental Data ####
#bifurcation index is from Mike Malick and can be downloaded from the share directory on his 
#Github here: https://github.com/michaelmalick/bifurcation-index 
#see: Malick, M.J., et al. 2017. 
#Effects of the North Pacific Current on the productivity of 163 Pacific salmon stocks. 
#Fisheries Oceanography 26:268--281. https://doi.org/10.1111/fog.12190

Bifurcation<-read.csv('data/Environment/BifurcationIndex.csv')%>%
mutate(stand_BI=scale(BI))%>%
  rename(calendar_year=Year)

# We used phonological indices to summarize BEUTI at annual scales at 1 degree latitude scales
#data was provided by Ellen M. Jorgensen and was generated in:
# Jorgensen, E. M., Hazen, E. L., Jacox M. G., Buil, M. P., Schroeder, I., Bograd, S. J.
# Physical and biogeochemical phenology of the coastal upwelling in the California Current Syste
# Geophysical Research Letters. In review as of 06/17/2024

Beuti_tumi<-read.csv('data/Environment/BEUTI_Phen/beuti_tumi.csv')%>%
  dplyr::select(year, X47N, X45N, X42N,X39N)%>%
  pivot_longer(c(X47N,X45N,X42N,X39N),names_to = "location1",values_to="Beuti_tumi" )%>%
  mutate(location = ifelse(location1 == 'X47N', "47N",
                    ifelse(location1 == 'X45N', "45N",
                    ifelse(location1 == 'X42N',"42N",
                    ifelse(location1 == 'X39N', "39N",
                         location1)))))%>%
  dplyr::select(-location1)%>%
  group_by(location)%>%
  mutate(Beuti_tumi_scale=scale(Beuti_tumi))%>%
  rename(calendar_year=year)%>%
  ungroup()

Beuti_sti<-read.csv('data/Environment/BEUTI_Phen/beuti_sti.csv')%>%
  dplyr::select(year, X47N, X45N, X42N,X39N)%>%
  pivot_longer(c(X47N,X45N,X42N,X39N),names_to = "location1",values_to="Beuti_sti" )%>%
  mutate(location = ifelse(location1 == 'X47N', "47N",
                    ifelse(location1 == 'X45N', "45N",
                    ifelse(location1 == 'X42N',"42N",
                    ifelse(location1 == 'X39N', "39N",
                         location1)))))%>%
  dplyr::select(-location1)%>%
  group_by(location)%>%
  mutate(Beuti_sti_scale=scale(Beuti_sti))%>%
  rename(calendar_year=year)%>%
  ungroup()


# SST Data can be accessed directly from the NOAA ERDAPP server
# polygons were provided by Jennifer Gosselin and Caitlin O'Brien
# see 01_SSTdata.R for code to aggregate and compile SST data that produces
# the 'ERSST_poly.rds' file

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

SST_winter<-SST_seasonal%>%filter(season=='Winter')%>%rename(sst_winter=seasonal_sst)%>%dplyr::select(-season)
SST_spring<-SST_seasonal%>%filter(season=='Spring')%>%rename(sst_spring=seasonal_sst)%>%dplyr::select(-season)

#### basin scael indices ####
dataInfo <- rerddap::info('PDO')
npgo<- read.csv("data/Environment/OC_NPGO.csv")
NPGO<-npgo%>%
  add_column('year'=as.numeric(format(as.Date(npgo$time),"%Y")))%>%
  add_column('month'=as.numeric(format(as.Date(npgo$time),"%m")))%>%
  add_column(variable="NPGO", dataset="NPGO", location="Basin Scale")%>%
  mutate()

oni<- read.csv("data/Environment/OC_ONI.csv")
ONI<-oni%>%
  add_column('year'=as.numeric(format(as.Date(oni$time),"%Y")))%>%
  add_column('month'=as.numeric(format(as.Date(oni$time),"%m")))%>%
  add_column(variable="ONI", dataset="ONI", location="Basin Scale")%>%
  dplyr::select(!time)

pdo<- read.csv("data/Environment/PDO.csv")%>%
  filter(year>1970)
PDO<-pdo%>%
  mutate(pdo_winter=scale(rowMeans(pdo%>%dplyr::select('X1', 'X2', 'X3'))))%>%
  mutate(pdo_spring=scale(rowMeans(pdo%>%dplyr::select('X4', 'X5', 'X6'))))%>%
 dplyr::select(year, pdo_spring, pdo_winter)%>%
  rename(calendar_year=year)


#### Combing Environmental with Chinook####

Atlas<-data.frame(area=na.omit(unique(Marine_Survival$area)),
     region=c("Salish Sea", "Salish Sea", "Salish Sea", "Columbia-Snake", 
              "Oregon Coast","Salish Sea", "Salish Sea", "Salish Sea", 
              "Columbia-Snake", "Washington Coast", "Columbia-Snake","Salish Sea",
              "Salish Sea", "Salish Sea", "Oregon Coast","Rogue-Klamath"))

Survival_combined<-Marine_Survival%>%dplyr::left_join(SST_winter,by=c('ecoregion', 'calendar_year'))%>%
                                dplyr::left_join(SST_spring, by=c('ecoregion', 'calendar_year'))%>%
  dplyr::left_join(Beuti_sti, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Beuti_tumi, by=c('location', 'calendar_year'))%>%
  dplyr::left_join(Bifurcation, by=c('calendar_year'))%>%
    dplyr::left_join(PDO, by=c('calendar_year'))%>%
  dplyr::left_join(Atlas,by=c('area'))
write.csv(Survival_combined, file = "data/Salmon/Survival_combined.csv")

dd<-Survival_combined%>%
  dplyr::select(calendar_year,ecoregion,location,Beuti_sti_scale,Beuti_tumi_scale,stand_BI,sst_spring,sst_winter)%>%
  distinct()
##### Exploratory Plots #####
data<-read.csv("data/Salmon/Survival_combined.csv")
unique(data$stock_code)
TS.survival<-ggplot(data =Survival_combined,
              aes(x =calendar_year, y = log(Marine.Survival), group = stock_name, col=stock_name))+
  #group=region)) +
  facet_wrap(.~region, ncol = 2, labeller = label_wrap_gen(25), scales="free_y") +
  geom_line(aes(group = stock_name,col=stock_name))+
  scale_y_continuous(name ="Log Age-2 Survival" )+
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
pairs(~ stand_BI+
        sst_spring+
        Beuti_sti_scale+
        Beuti_tumi_scale+
        sst_winter,
      data = na.omit(dd), upper.panel = panel.cor,         # Disabling the upper panel
#      data = na.omit(dd%>%filter(ecoregion==2, location=='45N')), upper.panel = panel.cor,         # Disabling the upper panel
      diag.panel = panel.hist,
      lower.panel = panel.smooth)

dev.off()
vif.covs<-car::vif(lm(Marine.Survival~stand_BI+sst_spring+Beuti_sti_scale+Beuti_tumi_scale+
        sst_winter, data=Survival_combined))
write.table(vif.covs, file = "Results/Tables/VIF.txt", sep = ",", quote = FALSE, row.names = F)
