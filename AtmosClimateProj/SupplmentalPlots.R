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
library(plotrix)
#### Reading in data/assigning colors #####
col2<-pnw_palette("Sunset2",4,type="discrete")
col<-pnw_palette("Sunset2",3,type="discrete")
col3<-pnw_palette("Sunset2",8,type="continuous")

col<-pnw_palette("Sunset2",3,type="discrete")
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
climate_dat_cop <-readRDS(here('data/physical/climate_dat_cop.rds'))
bakunsites <- read.csv(here('data/physical/Bakun/MapLocations.csv'))%>%
  mutate(longitude=longitude+360)
sites <- st_as_sf(data.frame(bakunsites[,1:2]), coords = c("longitude","latitude"), crs = 4326, 
                  agr = "constant")

#### Making the Map #####
map<-ggplot() +
  geom_polygon(aes(x=c(-105+360, -113+360, -127+360,-105+360,-105+360),
                   y=c(22.1, 22.1,34.4486,34.4486,20)),
               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110+360, -127+360,-130+360,-110+360,-110+360), 
                   y=c(34.4486,34.4486,40.4401,40.4401,34.4486)),
               fill='white', col='black',alpha=0.6)+
  geom_polygon(aes(x=c(-110+360, -130+360,-130+360,-110+360,-110+360), 
                   y=c(40.4401,40.4401,49.5,49.5,40.4401)),
               fill='white', col='black',alpha=0.6)+
  geom_sf(data = world)+
  annotate("rect", xmin= -121.5+360, xmax = -109+360, ymin = 42, ymax = 48.8, 
           fill = 'white', col='black',size = 0.8, lwd=0.2) +
  geom_sf(fill='grey95') +
  geom_sf(data = sites, size = c(rep(2,68+35+16), rep(3,2)), 
          shape = c(rep(24,68), rep(21,35),rep(23,16),rep(22,2)), 
          col = c(rep('black',68+35+18)), 
          fill = c(rep(col3[3],68), rep(col3[7],35),rep(col3[5],16),rep(col3[8],2))) +
  coord_sf(xlim = c(-132+360, -108+360), ylim = c(22, 50), expand = FALSE)+
  ylab(" ")+
  xlab(" ")+
  annotation_scale()+
  annotation_north_arrow(which_north = "true",pad_x = unit(0.25, "in"), 
                                 pad_y = unit(0.25, "in"))+
  annotate(geom = "text", x = c(-127+360,-118.5+360,-129+360), y = c(37,28,45), 
           label = str_wrap(c("Central", "Southern","Northern"), width = 20),
           fontface = "italic", color = "grey22", size = 3.75, angle=c('290', '311','270')) +
  annotate(geom = "text", x = c(-114+360,-114+360,-114+360,-114+360,-114+360), y = c(48,46.5,45, 44,43), 
           label = str_wrap(c("Bakun Index","CC Regions", "Newport Line","CalCOFI", "RREAS"), width = 22),
           color = "grey22", size =3.5) +
  annotate(geom = "text", x = c(-121+360,-117+360), y = c(41,35), 
           label = str_wrap(c("Cape Mendocino","Point Conception"), width = 20),
           fontface = "italic", color = "grey22", size = 3) +
  annotate("rect", xmin= -121+360, xmax = -119+360, ymin = 46, ymax = 47, 
           fill = 'white', col='black',size = 0.8, lwd=0.5) +
  annotate("line", x= c(-124.1+360, -124.65+360), y = c(44.652, 44.652),col=col2[1],size = 0.8, lwd=1) +
  annotate("line", x= c(-120.5+360, -119.5+360), y = c(45, 45),col=col2[1],size = 0.8, lwd=1) +
  
   theme(panel.background = element_rect(fill = "lightsteelblue2"),
        panel.border = element_rect(fill = NA),panel.grid.major = element_line(colour = "transparent"))

 map 
 pdf(file = "Output/Figures/Map.pdf",width = 6,height = 6)
 map
 dev.off()  

 
 ##### Reding in Bakun data ####
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
  add_column("Day"=as.numeric(format(as.Date(bakundat$time),"%d")))%>%
  add_column("YearDay"=as.numeric(yday(format(as.Date(bakundat$time)))))

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
bakun_daily <- bakun_region%>%
  group_by( station_id,Year, YearDay,region)%>%
  summarise(upwelling_index_sum = mean(na.omit(upwelling_index)))


##### Generating Cumulative Upwelling Plots #####

bakun_cum <- bakun_daily%>%
  group_by(station_id,Year,region)%>%
  reframe(upwelling_index_cum = cumsum(upwelling_index_sum))%>%
  add_column(YearDay=bakun_daily$YearDay)
period2=data.frame(period2=c('1967 - 1988', '1989 - 2013', '2014 - 2022'), period=c('1','2','3'))

bakun_time <-bakun_cum%>%
  filter(Year>1963 & Year<1989)%>%
  mutate(period='1')%>%
  bind_rows(bakun_cum%>%
              filter(Year>1989 & Year<2014)%>%
              mutate(period='2'))%>%
  bind_rows(bakun_cum%>%
              filter(Year>2013)%>%
              mutate(period='3'))%>%
  left_join(period2)


ggplot(data = bakun_time, aes(x = YearDay, y = upwelling_index_cum, col=Year)) +
  facet_wrap(.~station_id, ncol = 3, scales='free') +
  geom_line(size=0.75,aes(group=as.numeric(Year))) +
  # geom_smooth(col='red',aes(group=period)) +
  #scale_x_continuous(name = "Day") +
  #scale_color_manual(values=col3)+
  theme_bw()

bakun_time[which.min(bakun_time$upwelling_index_cum),]

CUM<-ggplot(data = bakun_time%>%filter(station_id!='60N'), aes(x = YearDay, y = upwelling_index_cum)) +
  facet_wrap(station_id~region, scales='free', ncol=3) +
  #  ggtitle(paste(bakun_time$station_id, "/n", bakun_time$region))+
  geom_line(size=0.75,aes(group=as.numeric(Year)),col='grey') +
  geom_smooth(aes(group=period2, colour=period2)) +
  #scale_x_continuous(name = "Day") +
  scale_color_manual(values=col[1:3], name="Period")+
  xlab("Julian Day")+
  ylab(expression('CUI ' ~ m^3 ~ '/ s / 100m'))+
  theme_bw()

pdf(file = "Output/Figures/CumulativeUpwelling.pdf",   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 8)
CUM
#thetaplot
dev.off()

#### Looking at Correlation Across Sites #####
cov.dat.test <- bakun%>%select(station_id, Month, Day, Year, upwelling_index)%>%
  group_by(Month, Day, Year, station_id) %>%
  summarise(upwelling_index = mean(upwelling_index))%>%
  pivot_wider(names_from = "station_id", values_from = "upwelling_index")%>%
  ungroup()
cov.dat.test <- cov.dat.test[complete.cases(cov.dat.test),]%>%
  select(`24N`,`27N`,`30N`,`33N`,`36N`,`39N`,`42N`,`45N`,`48N`,`51N`,`54N`,`57N`,`60N`)
cov_result<-cov(cov.dat.test)
cor_result<-cor(cov.dat.test)

color2D.matplot(cov_result, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                #vcex = 2,
                vcol = "black",
                extremes = c("red", "yellow"))
axis(3, at = seq_len(ncol(cov_result)) - 0.5,
     labels = colnames(cov_result), tick = FALSE)
axis(2, at = seq_len(nrow(cov_result)) -0.5,
     labels = rev(rownames(cov_result)), tick = FALSE, las = 1)

color2D.matplot(cor_result, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                #vcex = 2,
                vcol = "black",
                extremes = c("white","darkgreen"))
axis(3, at = seq_len(ncol(cov_result)) - 0.5,
     labels = colnames(cov_result), tick = FALSE)
axis(2, at = seq_len(nrow(cov_result)) -0.5,
     labels = rev(rownames(cov_result)), tick = FALSE, las = 1)
axis(1, at = c(2,5,7.5,10.5),
     labels = c("SCC","CCC", "NCC", "GoA"), tick = FALSE, las = 1)
axis(1, at = c(0,4,6,9,12),lwd = 4,labels=FALSE, tick = TRUE, las = 1)
axis(4, at = c(10.5,7.5,5,2)-0.5,
     labels = c("SCC","CCC", "NCC", "GoA"), tick = FALSE, las = 2,hadj=.5)
axis(4, at = c(0,3,6,8,12),lwd = 4,labels=FALSE, tick = TRUE, las = 2)
abline(h=c(0,3,6,8,12),lwd = 3)
abline(v=c(0,4,6,9,12),lwd = 3)
segments(x0=c(0,4,6,9,12),y0=c(0,0,0,0,0),x1=c(0,4,6,9,12),y1=c(12,12,12,12,12),col="black", lwd=3)
segments(x0=c(0,0,0,0,0),y0=c(0,3,6,8,12),x1=c(12,12,12,12,12),y1=c(0,3,6,8,12),col="black", lwd=3)


#### Looking at data by station and Month ####
monthly <- readRDS(here("data/physical/climate_dat_monthly.rds"))

bakun_monthly <- bakun_region%>%
  group_by( station_id,Year, Month,region)%>%
  summarise(upwelling_index_sum = mean(na.omit(upwelling_index)))%>%
  merge(monthly, by=c('Year', 'Month', 'region'))

reg.names=unique(monthly$region)
stat.names=unique(bakun_monthly$station_id)
figure<-function(dat){
  ggplot(data = dat, aes(y = upwelling_index_sum, x = PDO,col=period)) +
    facet_wrap(.~Month, ncol = 4, scales='free') +
    geom_point(size=0.75,aes(col=period)) +
    ggtitle(dat$station_id[1])+
    geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
    scale_y_continuous(name = "Upwelling (Bakun 1˚ 6-hourly)") +
    scale_x_continuous(name = "PDO") +
    scale_color_manual(values =  col[1:3],name="Period", labels=c('1967 - 1988', '1989 - 2013', '2014 - 2022'))+
    theme(plot.title = element_text(hjust = 0.5))+
    theme_bw()
}

figure(bakun_monthly%>%filter(station_id=='36N'))
figure(bakun_monthly%>%filter(station_id=='39N'))
figure(bakun_monthly%>%filter(station_id=='42N'))
figure(bakun_monthly%>%filter(station_id=='45N'))
figure(bakun_monthly%>%filter(station_id=='48N'))


figure(bakun_monthly%>%filter(region=="Northern CC"))
for(i in 1:length(stat.names)){
  print(figure(bakun_monthly%>%filter(station_id==stat.names[i])))
  #print(figure(df.levels%>%filter(Covar.Name==names.covars[i])))
}

pdf(file = "Output/Figures/Station-Month/36N.pdf",width = 6,height = 6)
figure(bakun_monthly%>%filter(station_id=='36N'))
dev.off()
pdf(file = "Output/Figures/Station-Month/39N.pdf",width = 6,height = 6)
figure(bakun_monthly%>%filter(station_id=='39N'))
dev.off()
pdf(file = "Output/Figures/Station-Month/42N.pdf",width = 6,height = 6)
figure(bakun_monthly%>%filter(station_id=='42N'))
dev.off()
pdf(file = "Output/Figures/Station-Month/45N.pdf",width = 6,height = 6)
figure(bakun_monthly%>%filter(station_id=='45N'))
dev.off()
pdf(file = "Output/Figures/Station-Month/48N.pdf",width = 6,height = 6)
figure(bakun_monthly%>%filter(station_id=='48N'))
dev.off()
pdf(file = "Output/Figures/Station-Month/ALLStations.pdf",width = 6,height = 6)
for(i in 1:length(stat.names)){
  print(figure(bakun_monthly%>%filter(station_id==stat.names[i])))
  #print(figure(df.levels%>%filter(Covar.Name==names.covars[i])))
}
dev.off()



#### SLP Correlation Analysis ####
correlation <-readRDS(here('data/physical/correlation_analysis.rds'))
correlation2 <-readRDS(here('data/physical/correlation_analysis_indices.rds'))
correlationdiff <-readRDS(here('data/physical/correlation_analysis_diff.rds'))


col<-pnw_palette("Sunset2",3,type="discrete")
slp_cor <- correlation2%>%
  filter(var=="SLP")


slp.plot <-ggplot() + 
  geom_raster(data=slp_cor, aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_wrap(~analysis, ncol = 3) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  scale_fill_gradient2(low = "blue", high = "red", name = "Coefficient \n (slope)") + 
  ggtitle("Spring SLP Anomalies (Pa) vs. Climate Indices")+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks = c(190, 210, 230))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))

pdf(file = "Output/Figures/SLP.pdf",width = 7,height = 7)
slp.plot
dev.off()

#### SST Correlation Analysis ####

col<-pnw_palette("Sunset2",3,type="discrete")
sst_cor <- correlation2%>%
  filter(var=="SST")


sst.plot <-ggplot() + 
  geom_raster(data=sst_cor, aes(x=longitude,y=latitude,fill = coefficient)) + 
  facet_wrap(~analysis, ncol = 3) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(190,240), ylim=c(30,60)) +
  scale_fill_gradient2(low = "blue", high = "red", name = "Coefficient \n (slope)") + 
  ggtitle("Spring SST Anomalies vs. Climate Indices")+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks = c(190, 210, 230))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))

pdf(file = "Output/Figures/SST.pdf",width = 7,height = 7)
sst.plot
dev.off()
