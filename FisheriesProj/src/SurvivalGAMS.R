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

Survival_combined<- read.csv("data/Survival_combined.csv")

##### Fitting Temporal GAMS #####
dd <- Survival_combined[complete.cases(Survival_combined), ]
dd <-cbind(dd,area1=as.factor(dd$area))

#smoothed term by year with random effect smoothed term
t1<-gam(Marine.Survival ~ s(calendar_year, k=6, m=2)+
              s(calendar_year, area1, k=6, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
draw(t1)

#smoothed term by year for each region (unshared)
t2<-gam(Marine.Survival ~s(calendar_year, by=area1,k=6),
data=dd, method="REML")
#plot.gam(g1)
draw(t2)

##Single smoothed term by year
t3<-gam(Marine.Survival ~ s(calendar_year, k=6),
data=dd, method="REML")
#plot.gam(g1)
draw(t3)

#shared smooth term by year with intercept random effect by region but common smoothed term
t4<-gam(Marine.Survival ~ s(calendar_year, k=5, bs="tp") +
s(area1, k=12, bs="re"),
data=dd, method="REML", family="gaussian")
draw(t4)

tAIC<-c(AIC(t1),AIC(t2),AIC(t3),AIC(t4))

##### Fitting Covariate GAMS #####
g1<-gam(Marine.Survival ~ s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
#plot.gam(g1)
draw(g1)

g2<-gam(Marine.Survival ~ s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g3)
draw(g2)

g3<-gam(Marine.Survival ~ s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g3)
draw(g3)

g2.3<-gam(Marine.Survival ~ s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g2.3)

g1.2<-gam(Marine.Survival ~ s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g1.2)

g1.3<-gam(Marine.Survival ~ s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g1.3)

g1.2.3<-gam(Marine.Survival ~s(Beuti_spring, k=5, m=2)+
              s(Beuti_spring, area1, k=5, bs="fs", m=2)+
              s(seasonal_STI, k=5, m=2)+
              s(seasonal_STI, area1, k=5, bs="fs", m=2)+
            s(sst_spring, k=5, m=2)+
              s(sst_spring, area1, k=5, bs="fs", m=2),
data=dd, method="REML")
summary(g1.2.3)
draw(g1.2)

#### Creating AIC Table & Predictions ####
min=AIC(g1.2)
aic <- data.frame(cbind(Model=c("Spring SST", "STI", "Spring BEUTI", "Spring SST; Spring BEUTI",
                               "STI; Spring BEUTI", "Spring SST; STI",
                     "Spring SST; Spring BEUTI; STI"),
  AIC=c(AIC(g1),AIC(g2),AIC(g3), AIC(g1.2),AIC(g1.3), AIC(g2.3), AIC(g1.2.3)),
  delAIC=c(AIC(g1)-min,AIC(g2)-min,AIC(g3)-min, AIC(g1.2)-min,AIC(g1.3)-min, AIC(g2.3)-min, AIC(g1.2.3)-min)))
aic
write.table(aic, file = "Results/Tables/AIC.txt", sep = ",", quote = FALSE, row.names = F)
# setup prediction data
Beuti_pred <- with(dd,
                      expand.grid(Beuti_spring=seq(min(Beuti_spring), max(Beuti_spring), length=50),
                                  sst_spring=rep(0, length=50),
                                  #stand_BI=seq(min(stand_BI), max(stand_BI), length=50),
                                  area1=levels(area1)))


SST_pred <- with(dd,
                      expand.grid(sst_spring=seq(min(sst_spring), max(sst_spring), length=50),
                                  Beuti_spring=rep(0, length=50),
                                  #stand_BI=seq(min(stand_BI), max(stand_BI), length=50),
                                  area1=levels(area1)))


#### Generating the plots ####
p1 <- draw(g1.2, select = "s(Beuti_spring)")+
  xlab("Spring BEUTI")+
  ggtitle("Fixed Effects") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p1
p2 <- draw(g1.2, select = "s(Beuti_spring,area1)")+
  xlab("Spring BEUTI")+
  ggtitle("Random Effects") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+ 
  guides(color =guide_legend(title="RMIS Region"))
p2
p3 <- draw(g1.2, select = "s(sst_spring)")+
  xlab("Spring SST")+
  ggtitle("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
p3
p4 <- draw(g1.2, select = "s(sst_spring,area1)")+
  xlab("Spring SST")+
  ggtitle("") +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5),
    legend.text = element_text(color = "white"),
    legend.title = element_text(color = "white"),
    legend.key = element_rect(fill = "white")) + 
  guides(color = guide_legend(override.aes= list(alpha = 0, color = "white"))) +
  theme(legend.key=element_rect(colour="white"))
p4
p1 + p2 + p3 + p4+plot_layout(ncol = 2,nrow=2)

# make the prediction, add this and a column of standard errors to the prediction
# data.frame. Predictions are on the log scale.
dd2<-cbind(dd$area1, dd$stock_code)
Beuti_pred<- cbind(Beuti_pred,
                       predict(g1.2, 
                               Beuti_pred, 
                               se.fit=TRUE, 
                               type="response"))
SST_pred<- cbind(SST_pred,
                       predict(g1.2, 
                               SST_pred, 
                               se.fit=TRUE, 
                               type="response"))
#g1.2_pred <-merge(g1.2_pred,dd2)
# make the plot. Note here the use of the exp() function to back-transform the
# predictions (which are for log-uptake) to the original scale
Beuti_pred_plot<-ggplot(data=dd, aes(y=Marine.Survival, x=Beuti_spring, group=area1)) +
  facet_wrap(~area1,scales="free_y") +
   geom_ribbon(aes(ymin=fit - 2*se.fit, ymax=fit + 2*se.fit, x=Beuti_spring),
              data=Beuti_pred, 
              alpha=0.3, 
              inherit.aes=FALSE) +
  geom_point(aes(col=stock_code))+#
  geom_line(aes(y=fit), data=Beuti_pred)


SST_pred_plot<-ggplot(data=dd, aes(y=Marine.Survival, x=sst_spring, group=area1)) +
  facet_wrap(~area1,scales="free_y") +
   geom_ribbon(aes(ymin=fit - 2*se.fit, ymax=fit + 2*se.fit, x=sst_spring),
              data=SST_pred, 
              alpha=0.3, 
              inherit.aes=FALSE) +
  geom_point(aes(col=stock_code))+#
  geom_line(aes(y=fit), data=SST_pred)
  
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
