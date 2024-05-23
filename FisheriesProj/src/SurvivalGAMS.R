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
t1<-gam(Marine.Survival ~ s(calendar_year, k=6),
data=dd, method="REML")
#plot.gam(g1)
#draw(t1)

#smoothed term by year for each region (unshared)
t2<-gam(Marine.Survival ~s(calendar_year, by=area1,k=6),
data=dd, method="REML")
#plot.gam(g1)
#draw(t2)

#smoothed term by year with random effect smoothed term
t3<-gam(Marine.Survival ~ s(calendar_year, k=6, m=2)+
              s(calendar_year, area1, k=6, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
#draw(t3)

tAIC<-c(AIC(t1),AIC(t2),AIC(t3))
mintAIC=min(tAIC)
taic <- data.frame(cbind(Model=c("Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term"),
  AIC=c(AIC(t1),AIC(t2),AIC(t3)),
  delAIC=c(AIC(t1)-mintAIC,AIC(t2)-mintAIC,AIC(t3)-mintAIC)))
taic
write.table(taic, file = "Results/Tables/tAIC.txt", sep = ",", quote = FALSE, row.names = F)
##### Fitting Beuti GAMS #####
dd <- Survival_combined[complete.cases(Survival_combined), ]
dd <-cbind(dd,area1=as.factor(dd$area))
##Single smoothed term by year
b1<-gam(Marine.Survival ~ s(Beuti_spring, k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(b1)

#smoothed term by year for each region (unshared)
b2<-gam(Marine.Survival ~s(Beuti_spring, by=area1,k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(b2)

#smoothed term by year with random effect smoothed term
b3<-gam(Marine.Survival ~ s(Beuti_spring, k=4, m=2)+
              s(Beuti_spring, area1, k=4, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
#draw(b3)

bAIC<-c(AIC(b1),AIC(b2),AIC(b3))
minbAIC=min(bAIC)
baic <- data.frame(cbind(Model=c("Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term"),
  AIC=c(AIC(b1),AIC(b2),AIC(b3)),
  delAIC=c(AIC(b1)-minbAIC,AIC(b2)-minbAIC,AIC(b3)-minbAIC)))
baic

###### Fitting SST GAMS #####

##Single smoothed term by year
s1<-gam(Marine.Survival ~ s(sst_spring, k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(s1)

#smoothed term by year for each region (unshared)
s2<-gam(Marine.Survival ~s(sst_spring, by=area1,k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(s2)

#smoothed term by year with random effect smoothed term
s3<-gam(Marine.Survival ~ s(sst_spring, k=4, m=2)+
              s(sst_spring, area1, k=4, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
#draw(s3)

sAIC<-c(AIC(s1),AIC(s2),AIC(s3))
minsAIC=min(sAIC)
saic <- data.frame(cbind(Model=c("Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term"),
  AIC=c(AIC(s1),AIC(s2),AIC(s3)),
  delAIC=c(AIC(s1)-minsAIC,AIC(s2)-minsAIC,AIC(s3)-minsAIC)))
saic

##### Fitting STI GAMS #####

##Single smoothed term by year
st1<-gam(Marine.Survival ~ s(seasonal_STI, k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(st1)

#smoothed term by year for each region (unshared)
st2<-gam(Marine.Survival ~s(seasonal_STI, by=area1,k=4),
data=dd, method="REML")
#plot.gam(g1)
#draw(st2)

#smoothed term by year with random effect smoothed term
st3<-gam(Marine.Survival ~ s(seasonal_STI, k=4, m=2)+
              s(seasonal_STI, area1, k=4, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
#draw(st3)

stAIC<-c(AIC(st1),AIC(st2),AIC(st3))
minstAIC=min(stAIC)
staic <- data.frame(cbind(Model=c("Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term"),
  AIC=c(AIC(st1),AIC(st2),AIC(st3)),
  delAIC=c(AIC(st1)-minstAIC,AIC(st2)-minstAIC,AIC(st3)-minstAIC)))
staic


##### Fitting Covariate GAMS #####
g.st.s<-gam(Marine.Survival ~ s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g.st.s)

gb.s<-gam(Marine.Survival ~ s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(gb.s)

gst.b<-gam(Marine.Survival ~ s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(gst.b)

gst.b.s<-gam(Marine.Survival ~s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2)+
              s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(gst.b.s)
#draw(gst.b.s)

#### testing return season ####

g.ra<-gam(Marine.Survival ~ s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2)+
              s(Beuti_spring,season1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring,area1, k=5, bs="fs", m=2)+
              s(sst_spring,season1, k=5, bs="fs", m=2),
data=dd, method="REML")
draw(g.ra)
summary(g.ra)

#### Creating AIC Table & Predictions ####
gAIC=c(AIC(b1),AIC(b2),AIC(b3),
       AIC(s1),AIC(s2),AIC(s3),
       AIC(st1),AIC(st2),AIC(st3),
       AIC(g.st.s),AIC(gb.s), AIC(gst.b), AIC(gst.b.s), AIC(g.ra))
min=min(gAIC)
aic <- data.frame(cbind(Model=c("Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term",
                     "Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term",
                     "Single Shared Smooth Term", "Unshared Smooth Term",
                     "Single Shared Smooth Term and Group-Level Smooth Term",
                      "Single Shared Smooth Term and Group-Level Smooth Term",
                      "Single Shared Smooth Term and Group-Level Smooth Term",
                       "Single Shared Smooth Term and Group-Level Smooth Term",
                        "Single Shared Smooth Term and Group-Level Smooth Term",
                        "Single Shared Smooth Term and Group-Level (RMIS & Season) Smooth Term"),
                     Covariates=c("BEUTI", "BEUTI","BEUTI",
                                  "SST","SST","SST",
                                  "STI","STI","STI",
                                  "SST & STI", "BEUTI & SST", "STI & BEUTI", 
                                  "STI & BEUTI & SST", "STI & SST"),
  AIC=gAIC,
  delAIC=gAIC-min))
aic
write.table(aic, file = "Results/Tables/AIC.txt", sep = ",", quote = FALSE, row.names = F)
# setup prediction data
Beuti_pred <- with(dd,
                      expand.grid(Beuti_spring=seq(min(Beuti_spring), max(Beuti_spring), length=50),
                                  sst_spring=rep(0, length=50),
                                  #stand_BI=seq(min(stand_BI), max(stand_BI), length=50),
                                  area1=levels(area1),season1=levels(season1)))


SST_pred <- with(dd,
                      expand.grid(sst_spring=seq(min(sst_spring), max(sst_spring), length=50),
                                  Beuti_spring=rep(0, length=50),
                                  #stand_BI=seq(min(stand_BI), max(stand_BI), length=50),
                                  area1=levels(area1),season1=levels(season1)))


#### Generating the plots ####
p1 <- draw(g.ra, select = "s(Beuti_spring)")+
  xlab("Beuti_spring")+
  ggtitle("Global Smooth") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
p2 <- draw(g.ra, select = "s(Beuti_spring,area1)")+
  xlab("Beuti_spring")+
  ggtitle("RMIS Group-Level") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  guides(color =guide_legend(title="RMIS Region"))
p2
p3 <- draw(g.ra, select = "s(Beuti_spring,season1)")+
  xlab("Beuti_spring")+
  ggtitle("Run Type Group-Level") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  guides(color =guide_legend(title="Run Type"))
p3
p4 <- draw(g.ra, select = "s(sst_spring)")+
  xlab("Spring SST")+
  ggtitle("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p4
p5 <- draw(g.ra, select = "s(sst_spring,area1)")+
  xlab("Spring SST")+
  ggtitle("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  guides(color = guide_legend(override.aes= list(alpha = 0, color = "white"))) +
  theme(legend.key=element_rect(colour="white"))
p5
p6 <- draw(g.ra, select = "s(sst_spring,season1)")+
  xlab("Spring SST")+
  ggtitle("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  guides(color = guide_legend(override.aes= list(alpha = 0, color = "white"))) +
  theme(legend.key=element_rect(colour="white"))
p6
p1 + p2 + p3 + p4+ p5 + p6+plot_layout(ncol = 3,nrow=2)

# make the prediction, add this and a column of standard errors to the prediction
# data.frame. Predictions are on the log scale.
dd2<-cbind(dd$area1,dd$season1, dd$stock_code)
Beuti_pred<- cbind(Beuti_pred,
                       predict(g.ra, 
                              Beuti_pred, 
                               se.fit=TRUE, 
                               type="response"))
SST_pred<- cbind(SST_pred,
                       predict(g.ra, 
                               SST_pred, 
                               se.fit=TRUE, 
                               type="response"))
pred.b<-unique(dd%>%dplyr::select(area1, season1))%>%left_join(Beuti_pred)
pred.s<-unique(dd%>%dplyr::select(area1, season1))%>%left_join(SST_pred)
p<-unique(dd%>%dplyr::select(stock_code, Code2))

#g1.2_pred <-merge(g1.2_pred,dd2)
# make the plot. Note here the use of the exp() function to back-transform the
# predictions (which are for log-uptake) to the original scale
Beuti_pred_plot<-ggplot(data=dd, aes(y=Marine.Survival, x=Beuti_spring, group=area1)) +
  facet_wrap(~area1,scales="free_y") +
   geom_ribbon(aes(ymin=fit - 2*se.fit, ymax=fit + 2*se.fit, 
                   x=Beuti_spring,group=season1, fill=season1),
              data=pred.b, 
              alpha=0.2, 
              inherit.aes=FALSE) +
  geom_point(aes(col=season1, group=season1,shape=stock_code))+#
  scale_shape_manual(values=p$Code2)+
  geom_line(aes(y=fit, group=season1, col=season1), data=pred.b)+
  guides(shape = FALSE)
Beuti_pred_plot

SST_pred_plot<-ggplot(data=dd, aes(y=Marine.Survival, x=sst_spring, group=area1)) +
  facet_wrap(~area1,scales="free_y") +
   geom_ribbon(aes(ymin=fit - 2*se.fit, ymax=fit + 2*se.fit, 
                   x=sst_spring,group=season1, fill=season1),
              data=pred.s, 
              alpha=0.1, 
              inherit.aes=FALSE) +
    geom_point(aes(col=season1, group=season1,shape=stock_code))+#
  scale_shape_manual(values=p$Code2)+
  geom_line(aes(y=fit, group=season1, col=season1), data=pred.s)+
  guides(shape = FALSE)
  SST_pred_plot
pdf(file = "Results/Figures/SST_pred_plot.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 7)
SST_pred_plot
dev.off()

pdf(file = "Results/Figures/Beuti_pred_plot.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 7)
Beuti_pred_plot
dev.off()

pdf(file = "Results/Figures/BestModel_plot.pdf",   # The directory you want to save the file in
    width = 11, # The width of the plot in inches
    height = 8.5)
p1 + p3 + p2 + p4 + p6+ p5+plot_layout(ncol = 3,nrow=2)

dev.off()




df_wider <- fem4mat %>%
  pivot_wider(names_from = brood_year, values_from = m, values_fill = list(m = NA)) %>%
  as.data.frame()

numeric_matrix <- as.matrix(df_wider[, -1])
