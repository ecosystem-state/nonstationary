
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


schroeder.nph <-read.csv('data/physical/year_mon_area_max_x_y_new.csv')%>%
  mutate(era= ifelse(Year<=1988,1, ifelse(Year>2012,3,2)))%>%
  mutate(Year_win = if_else(Month == 11|Month ==12, Year+1, Year))%>%
  mutate(dec.yr =as.numeric(as.character(Year)) + (as.numeric(Month)-0.5)/12)%>%
  mutate(area.anom = Area-mean(Area))
plot(schroeder.nph$dec.yr, schroeder.nph$Area, type="l")
plot(schroeder.nph$dec.yr, schroeder.nph$Max, type="l")
# set up as time series for breakpoint test
sst.dat <- ts(data=win.anom[2:70], 1951, 2019, frequency=1)




ggplot(data=schroeder.nph%>%filter(Month==4|Month==5|Month==6),
       aes(x,y, label = Year,col=as.factor(era)))+
  facet_wrap(.~Month, ncol = 3, scales='free') +
  geom_point()+
  geom_text()+
  stat_ellipse(level = 0.9) +
  theme_bw()
ggplot(data=schroeder.nph,
       aes(dec.yr,Max))+
  facet_wrap(.~Month, ncol = 3, scales='free') +
  geom_line()+
  theme_bw()


spring.schroeder <- schroeder.nph%>%
  filter(Month==4|Month==5|Month==6)%>%
  group_by(Year)%>%
  summarise(mean.x=mean(x), mean.y=mean(y),mean.max=mean(Max),mean.area=mean(Area),
            mean.area.anom=mean(area.anom))%>%
  mutate(era= ifelse(Year<=1988,1, ifelse(Year>2012,3,2)))%>%
  mutate(xy=sqrt(mean.x^2+mean.y^2), era=as.factor(era),
         era.lab = ifelse(era==1, '1967 - 1988', ifelse(era==2, "1989 - 2012", "2013 - 2023")))

means <- spring.schroeder%>%
  group_by(era.lab)%>%
  summarise(x=mean(mean.x),y=mean(mean.y),sd.x=sd(mean.x), sd.y=sd(mean.y), Year=0)%>%
  rename(mean.x=x, mean.y=y)
theme_set(theme_classic())

col<-pnw_palette("Sunset2",3,type="discrete")
a.plot <-ggplot(data=spring.schroeder,aes(mean.x,mean.y, label=Year,group=era.lab,col=era.lab))+
  geom_text(col='grey')+
  ggtitle("Center of North Pacific High") +
  geom_point(data=means)+
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  geom_errorbar(data=means,aes(ymin = mean.y-sd.y, ymax=  mean.y+sd.y), width=0.5) +
  geom_errorbar(data=means,aes(xmin = mean.x-sd.x, xmax=  mean.x+sd.x), width=0.5) +
  scale_colour_manual(values = c(col[1], col[2], col[3]), name = "") +
  ylab('y')+
  theme_bw() +
  xlab('x')
a.plot

plot(spring.schroeder$Year, spring.schroeder$xy, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.y, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.x, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.max, type="l")
plot(spring.schroeder$Year, spring.schroeder$mean.area.anom, type="l")

max.ts <- ts(data=spring.schroeder%>%select(mean.y), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max) 

max.ts <- ts(data=spring.schroeder%>%select(mean.max), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max) 

max.ts <- ts(data=spring.schroeder%>%select(mean.area.anom), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max) 

max.ts <- ts(data=spring.schroeder%>%select(xy), 1967, 2023, frequency=1)
# fit breakpoint model
bp.max <- breakpoints(max.ts ~ 1)
summary(bp.max) 


mod <- lm(mean.y ~ as.factor(era), data=spring.schroeder)

# and get predicted values of the model to plot
pred <- predict(mod, se=T, newdata = spring.schroeder)
spring.schroeder$mean <- pred$fit

# now save for a combined plot

# set the colors to use - colorblind pallette
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
theme_set(theme_classic())

b.plot <- ggplot(data=spring.schroeder, aes(Year, mean.y)) +
  geom_line(size=0.2) +
  xlab("")+
  geom_line(aes(Year, mean), color=cb[6], size=0.4) + 
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8,hjust = 0.5), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  ylab("Latitude") + ggtitle("Center of North Pacific High (y)") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  theme_bw() +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1966,2024)
b.plot


# and calculate standard deviation over 11-year rolling windows
NPH.max.sd <- rollapply(spring.schroeder%>%select(mean.max), 11, sd, fill=NA)
plot(1967:2023, NPH.max.sd , type="l")

# and calculate standard deviation over 11-year rolling windows
NPH.area.sd <- rollapply(spring.schroeder%>%select(mean.area.anom), 11, sd, fill=NA)
plot(1967:2023, NPH.area.sd, type="l")
# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.max.sd))
# fit the model
mod <- gam(mean.max ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

c.plot <- ggplot(plot.dat, aes(year, mean.max)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4)  +  theme_bw() +
  ylab("Standard deviation (pa)") +
  xlab("")+
  ggtitle("North Pacific High Intensity Variability") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
c.plot


# first, make a data frame
plot.dat <- data.frame(year=1972:2018, sd=na.omit(NPH.area.sd))
# fit the model
mod <- gam(mean.area.anom ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit  

d.plot <- ggplot(plot.dat, aes(year, mean.area.anom)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color=cb[6], size=0.4) +   theme_bw() +
  ylab(expression("Standard deviation "~(km^2))) +
  xlab("Year")+
  ggtitle("North Pacific High Area Variability") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  geom_vline(xintercept = 2012.5, lty=2, size=0.3) +
  xlim(1967,2020)
d.plot








png("figs/Fig 1.png", 6, 11, units="cm", res=300) 
ggarrange(a.plot, b.plot, c.plot, labels = c("A", "B", "C"),  
          font.label = list(size = 10, face="plain"), nrow=3, align="v", hjust = -1)
dev.off()

pdf("figs/Fig 1.pdf", 6/2.54, 11/2.54) 
ggarrange(a.plot, b.plot, c.plot, labels = c("A", "B", "C"),  
          font.label = list(size = 10, face="plain"), nrow=3, align="v", hjust = -1)
dev.off()