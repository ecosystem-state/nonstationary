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

Survival_combined<- read.csv("data/Salmon/Survival_combined.csv")%>%
  dplyr::select(-Mat.Rate,-season,-X)%>%
  dplyr::filter(stock_code!='IGHY'&stock_code!='TRHY')%>%
  mutate(area=as.factor(area), RunType=as.factor(RunType))%>%
  filter(ecoregion==1|ecoregion==2|ecoregion==3|ecoregion==4)%>%
  filter(calendar_year>=1988&calendar_year<=2023)
table(is.na(Survival_combined))["FALSE"]

saveRDS(Survival_combined, file = "data/Salmon/Survival_covariates.rds")

##### Fitting Temporal GAMS #####

##Single smoothed term by year
t1<-gam(log(Marine.Survival) ~ s(calendar_year, k=6),
        data=Survival_combined, method="REML")
plot.gam(t1)
draw(t1)

#smoothed term by year for each region (unshared)
t2<-gam(log(Marine.Survival) ~s(calendar_year, by=area,k=6),
        data=Survival_combined, method="REML")

draw(t2)

#smoothed term by year for each runtype (unshared)
t3<-gam(log(Marine.Survival) ~s(calendar_year, by=RunType,k=6),
        data=Survival_combined, method="REML")

draw(t3)

#smoothed term by year for each runtype (unshared) and region
t4<-gam(log(Marine.Survival) ~s(calendar_year, by=RunType,k=6)+
          s(calendar_year, by=area,k=6),
        data=Survival_combined, method="REML")

draw(t4)
#smoothed term by year with random effect smoothed term by runtype
t5<-gam(log(Marine.Survival) ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Survival_combined, method="REML")
#plot.gam(g1)
draw(t5)

#smoothed term by year with random effect smoothed term by Runtype
t6<-gam(log(Marine.Survival) ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2),
        data=Survival_combined, method="REML")
#plot.gam(g1)
draw(t6)


#smoothed term by year with random effect smoothed term
t7<-gam(log(Marine.Survival) ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Survival_combined, method="REML")
#plot.gam(g1)
draw(t7)


tAIC<-c(AIC(t1),AIC(t2),AIC(t3),AIC(t4),AIC(t5),AIC(t6),AIC(t7))
mintAIC=min(tAIC)
taic <- data.frame(cbind(Model=c("Single Shared Smooth Term", 
                                 "Unshared Smooth Term by Area",
                                 "Unshared Smooth Term by Runtype",
                                 "Unshared Smooth Term by Runtype and Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype and Area"
                                 ),
                         
                         AIC=tAIC,
                         delAIC=tAIC-min(tAIC)))
taic
#write.table(taic, file = "Results/Tables/tAIC.txt", sep = ",", quote = FALSE, row.names = F)


##### Time-varying variability ####

Survival_Variance<-Survival_combined%>%
  group_by(stock_code)%>%
 mutate(roll_sd=rollapply(log(Marine.Survival), 5, sd, fill=NA))%>%
  ungroup()

ggplot(aes(x=calendar_year,y=roll_sd),data=Survival_Variance)+
  geom_smooth()

##### Fitting Temporal GAMS #####

##Single smoothed term by year
V1<-gam(roll_sd ~ s(calendar_year, k=6),
        data=Survival_Variance, method="REML")

draw(V1)

#smoothed term by year for each region (unshared)
V2<-gam(roll_sd ~s(calendar_year, by=area,k=6),
        data=Survival_Variance, method="REML")

draw(V2)

#smoothed term by year for each runtype (unshared)
V3<-gam(roll_sd ~s(calendar_year, by=RunType,k=6),
        data=Survival_Variance, method="REML")

draw(V3)

#smoothed term by year for each runtype (unshared) and region
V4<-gam(roll_sd ~s(calendar_year, by=RunType,k=6)+
          s(calendar_year, by=area,k=6),
        data=Survival_Variance, method="REML")

draw(V4)
#smoothed term by year with random effect smoothed term by runtype
V5<-gam(roll_sd ~ s(calendar_year, k=5, m=2)+
          s(calendar_year, area, k=5, bs="fs", m=2),
        data=Survival_Variance, method="REML")
#plot.gam(g1)
draw(V5)

#smoothed term by year with random effect smoothed term by Runtype
V6<-gam(roll_sd ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2),
        data=Survival_Variance, method="REML")
#plot.gam(g1)
draw(V6)


#smoothed term by year with random effect smoothed term
V7<-gam(roll_sd ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Survival_Variance, method="REML")
#plot.gam(g1)
draw(V7)



VAIC<-c(AIC(V1),AIC(V2),AIC(V3),AIC(V4),AIC(V5),AIC(V6),AIC(V7))
Vaic <- data.frame(cbind(Model=c("Single Shared Smooth Term", 
                                 "Unshared Smooth Term by Area",
                                 "Unshared Smooth Term by Runtype",
                                 "Unshared Smooth Term by Runtype and Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype and Area"
),

AIC=VAIC,
delAIC=VAIC-min(VAIC)))
Vaic
write.table(taic, file = "Results/Tables/tAIC.txt", sep = ",", quote = FALSE, row.names = F)

##### Maturity through time #####

Mat_combined<- read.csv("data/Salmon/Survival_combined.csv")%>%
  dplyr::select(-season,-X)%>%
  dplyr::filter(stock_code!='IGHY'&stock_code!='TRHY')%>%
  mutate(area=as.factor(area), RunType=as.factor(RunType))%>%
  filter(ecoregion==1|ecoregion==2|ecoregion==3|ecoregion==4)%>%
  filter(calendar_year>=1988&calendar_year<=2023)
table(is.na(Mat_combined))["FALSE"]

##### Fitting Temporal GAMS #####

##Single smoothed term by year
t1<-gam(Mat.Rate ~ s(calendar_year, k=6),
        data=Mat_combined, method="REML")
draw(t1)

#smoothed term by year for each region (unshared)
t2<-gam(Mat.Rate ~s(calendar_year, by=area,k=6),
        data=Mat_combined, method="REML")

draw(t2)

#smoothed term by year for each runtype (unshared)
t3<-gam(Mat.Rate ~s(calendar_year, by=RunType,k=6),
        data=Mat_combined, method="REML")

draw(t3)

#smoothed term by year for each runtype (unshared) and region
t4<-gam(Mat.Rate~s(calendar_year, by=RunType,k=6)+
          s(calendar_year, by=area,k=6),
        data=Mat_combined, method="REML")

draw(t4)
#smoothed term by year with random effect smoothed term by runtype
t5<-gam(Mat.Rate ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Mat_combined, method="REML")
#plot.gam(g1)
draw(t5)

#smoothed term by year with random effect smoothed term by Runtype
t6<-gam(Mat.Rate~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2),
        data=Mat_combined, method="REML")
#plot.gam(g1)
draw(t6)


#smoothed term by year with random effect smoothed term
t7<-gam(Mat.Rate~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Mat_combined, method="REML")
draw(t7)


tAIC<-c(AIC(t1),AIC(t2),AIC(t3),AIC(t4),AIC(t5),AIC(t6),AIC(t7))
mintAIC=min(tAIC)
taic <- data.frame(cbind(Model=c("Single Shared Smooth Term", 
                                 "Unshared Smooth Term by Area",
                                 "Unshared Smooth Term by Runtype",
                                 "Unshared Smooth Term by Runtype and Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype and Area"
),

AIC=tAIC,
delAIC=tAIC-min(tAIC)))
taic
write.table(taic, file = "Results/Tables/tAIC.txt", sep = ",", quote = FALSE, row.names = F)


##### Time-varying variability ####

Mat_Variance<-Mat_combined%>%
  group_by(stock_code)%>%
  mutate(roll_sd=rollapply(Mat.Rate, 5, sd, fill=NA))%>%
  ungroup()

ggplot(aes(x=calendar_year,y=roll_sd),data=Survival_Variance)+
  geom_smooth()

##Single smoothed term by year
V1<-gam(roll_sd ~ s(calendar_year, k=6),
        data=Mat_Variance, method="REML")
draw(V1)

#smoothed term by year for each region (unshared)
V2<-gam(roll_sd ~s(calendar_year, by=area,k=6),
        data=Mat_Variance, method="REML")

draw(V2)

#smoothed term by year for each runtype (unshared)
V3<-gam(roll_sd ~s(calendar_year, by=RunType,k=6),
        data=Mat_Variance, method="REML")

draw(V3)

#smoothed term by year for each runtype (unshared) and region
V4<-gam(roll_sd~s(calendar_year, by=RunType,k=6)+
          s(calendar_year, by=area,k=6),
        data=Mat_Variance, method="REML")

draw(V4)
#smoothed term by year with random effect smoothed term by runtype
V5<-gam(roll_sd ~ s(calendar_year, k=6, m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Mat_Variance, method="REML")
#plot.gam(g1)
draw(V5)

#smoothed term by year with random effect smoothed term by Runtype
V6<-gam(roll_sd~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2),
        data=Mat_Variance, method="REML")
#plot.gam(g1)
draw(V6)


#smoothed term by year with random effect smoothed term
V7<-gam(roll_sd~ s(calendar_year, k=6, m=2)+
          s(calendar_year, RunType, k=6, bs="fs", m=2)+
          s(calendar_year, area, k=6, bs="fs", m=2),
        data=Mat_Variance, method="REML")
draw(V7)


VAIC<-c(AIC(V1),AIC(V2),AIC(V3),AIC(V4),AIC(V5),AIC(V6),AIC(V7))
Vaic <- data.frame(cbind(Model=c("Single Shared Smooth Term", 
                                 "Unshared Smooth Term by Area",
                                 "Unshared Smooth Term by Runtype",
                                 "Unshared Smooth Term by Runtype and Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Area",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype",
                                 "Single Shared Smooth Term and Group-Level Smooth Term by Runtype and Area"
),

AIC=VAIC,
delAIC=VAIC-min(VAIC)))
Vaic
