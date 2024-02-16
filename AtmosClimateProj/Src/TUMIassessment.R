library(ggrepel)
library(strucchange)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(mgcv)
library(strucchange)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# reading in the data 
TUMI <-read.csv('data/physical/cciea_OC_TUMI.csv')
TUMI <-TUMI %>%
  add_column('Year'=as.numeric(format(as.Date(TUMI$time),"%Y")))%>%
 #select(-time, -X48N)%>%
  pivot_longer(c(X33N,X36N, X39N,X42N,X45N,X48N), names_to = "Location", values_to = "TUMI")%>%
  mutate(Location=ifelse(Location=='X33N', '33N',
          ifelse(Location=='X36N','36N',
          ifelse(Location=='X39N', '39N',
          ifelse(Location=='X42N','42N',
          ifelse(Location=='X45N','45N',
          ifelse(Location=='X48N','48N',Location)))))))

p.TUMI<-ggplot(data=TUMI, 
       aes(x=Year,y=TUMI))+
  facet_wrap(.~Location, ncol = 3) +
  geom_line(col=cb[6])+
  geom_smooth(se=F, col='grey')+
  ggtitle('TUMI')+

  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
    geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  #geom_smooth(aes(group=era))+
  theme_bw()+
 theme(plot.title = element_text(hjust = 0.5))
  
y.ts <- ts(data=TUMI%>%filter(Location=="42N")%>%select(TUMI), 1967, 2022, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 

#### ROLLING WINDOW #####
TUMI.sd<-NA
Year<-unique(TUMI$Year)
loc<-unique(TUMI$Location)
for(i in 1:length(loc)){
  temp <- rollapply(TUMI%>%filter(Location==loc[i])%>%select(TUMI), 10, sd, fill=NA)
TUMI.sd<-cbind(TUMI.sd,temp)
}

TUMI.sd<-data.frame(TUMI.sd)
TUMI.sd<-TUMI.sd%>%
  select(-TUMI.sd)%>%
  mutate(Year=Year)
colnames(TUMI.sd)<-c(loc,"Year")
TUMI.sd<-TUMI.sd%>%
    pivot_longer(c("33N", "36N", "39N", "42N", "45N","48N"), names_to = "Location", values_to = "TUMI")

p.TUMI.sd<-ggplot(data=TUMI.sd, 
       aes(x=Year,y=TUMI))+
  facet_wrap(.~Location, ncol = 3) +
  geom_line(col=cb[6])+
  geom_smooth(se=F, col='grey')+
  ggtitle('TUMI SD 11-Year Rolling Window')+

  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
    geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  #geom_smooth(aes(group=era))+
  theme_bw()+
 theme(plot.title = element_text(hjust = 0.5))

TUMI.ts<-na.omit(TUMI.sd)  
y.ts <- ts(data=TUMI.ts%>%filter(Location=="42N")%>%select(TUMI), 1972, 2017, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 



# reading in the data 
STI <-read.csv('data/physical/cciea_OC_STI.csv')
STI <-STI %>%
  add_column('Year'=as.numeric(format(as.Date(STI$time),"%Y")))%>%
  select(-time)%>%
  pivot_longer(c(X33N,X36N, X39N,X42N,X45N), names_to = "Location", values_to = "STI")%>%
  mutate(Location=ifelse(Location=='X33N', '33N',
          ifelse(Location=='X36N','36N',
          ifelse(Location=='X39N', '39N',
          ifelse(Location=='X42N','42N',
          ifelse(Location=='X45N','45N',Location))))))

p.STI<-ggplot(data=STI, 
       aes(x=Year,y=STI))+
  facet_wrap(.~Location, ncol = 3) +
  geom_line(col=cb[6])+
  geom_smooth(se=F, col='grey')+
  ggtitle('STI')+
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
    geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  #geom_smooth(aes(group=era))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
  
y.ts <- ts(data=STI%>%filter(Location=="42N")%>%select(STI), 1967, 2022, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 


#### ROLLING WINDOW #####
STI.sd<-NA
Year<-unique(STI$Year)
loc<-unique(STI$Location)
for(i in 1:5){
  temp <- rollapply(STI%>%filter(Location==loc[i])%>%select(STI), 10, sd, fill=NA)
STI.sd<-cbind(STI.sd,temp)
}

STI.sd<-data.frame(STI.sd)
STI.sd<-STI.sd%>%
  select(-STI.sd)%>%
  mutate(Year=Year)
colnames(STI.sd)<-c(loc,"Year")
STI.sd<-STI.sd%>%
    pivot_longer(c("33N", "36N", "39N", "42N", "45N"), names_to = "Location", values_to = "STI")

p.STI.sd<-ggplot(data=STI.sd, 
       aes(x=Year,y=STI))+
  facet_wrap(.~Location, ncol = 3) +
  geom_line(col=cb[6])+
  geom_smooth(se=F, col='grey')+
  ggtitle('STI SD 11-Year Rolling Window')+

  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
    geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  #geom_smooth(aes(group=era))+
  theme_bw()+
 theme(plot.title = element_text(hjust = 0.5))

STI.ts<-na.omit(STI.sd)  
y.ts <- ts(data=STI.ts%>%filter(Location=="42N")%>%select(STI), 1972, 2017, frequency=1)
# fit breakpoint model
bp.y <- breakpoints(y.ts ~ 1)
summary(bp.y) 

### Printing plots
pdf("Output/Figures/Phenology/STI_RW.pdf", 7,5) 
p.STI.sd
dev.off()

pdf("Output/Figures/Phenology/STI.pdf", 7,5) 
p.STI
dev.off()

pdf("Output/Figures/Phenology/TUMI.pdf", 7,5) 
p.TUMI
dev.off()

pdf("Output/Figures/Phenology/TUMI_RW.pdf", 7,5) 
p.TUMI.sd
dev.off()


pdf("Output/Figures/Phenology/TUMI.pdf", 7,10) 
ggarrange( p.TUMI, p.TUMI.sd,labels = c("A", "B"),  
           font.label = list(size = 12, face="plain"), nrow=2)
dev.off()

pdf("Output/Figures/Phenology/STI.pdf", 7,10) 
ggarrange( p.STI, p.STI.sd,labels = c("A", "B"),  
           font.label = list(size = 12, face="plain"), nrow=2)
dev.off()
