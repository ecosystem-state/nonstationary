library(ncdf4)
library(chron)
library(tidyverse)
library(kohonen) # fitting
library(aweSOM) # plotting
library(SOMbrero) # plotting
library(paletteer) #colors
library(PNWColors) #more colors
library(here) #navigating folders
set.seed(1234)
here::i_am("som/sst_soms_spring.r")

#### Importing data ####

# start by loading NE Pacific SST/SLP - this may need to be moved into project if not locally stored
nc <- nc_open(here("som/copernicus_jun27.nc"))

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

#### SST Analysis ####

# get times and 
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$spring_year <- dates$year
dates$spring_year[which(dates$month %in% 7:12)] <- NA
dates$spring_year[which(dates$month %in% 1:3)] <- NA # incremenent 

var_name = "sst"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
X <- X[,-to_drop]# drop cells with NAs 
spring_lat <- lat[-to_drop]
spring_lon <- lon[-to_drop]
# Add block for spring average calculations
spring_years <- unique(dates$spring_year)
spring_years <- sort(spring_years[-which(is.na(spring_years))])
X_spring <- matrix(NA, length(spring_years), ncol(X))
for(i in 1:length(spring_years)) {
  X_spring[i,] = colMeans(X[which(dates$spring_year == spring_years[i]),])
}

#assigning names to time allows you to divide data by ERA
dimnames(X_spring) <- list(as.character(spring_years),paste("N", spring_lat, "E", spring_lon, sep=""))


climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))

CC.FMA <- climate_dat%>% filter(season=='Spring'&region=='Central CC')%>% 
  select(Year_lag, stand_bakun_seasonally) 
rownames(CC.FMA)<- CC.FMA$Year_lag

NCC.FMA <- climate_dat%>% filter(season=='Spring'&region=='Northern CC')%>% 
  select(Year_lag, stand_bakun_seasonally) 
rownames(NCC.FMA)<- NCC.FMA$Year_lag

SCC.FMA <- climate_dat%>% filter(season=='Spring'&region=='Southern CC')%>% 
  select(Year_lag, stand_bakun_seasonally) 
rownames(SCC.FMA)<- SCC.FMA$Year_lag

GoA.FMA <- climate_dat%>% filter(season=='Spring'&region=='GoA')%>% 
  select(Year_lag, stand_bakun_seasonally) 
rownames(GoA.FMA)<- GoA.FMA$Year_lag

SST1 <- X_spring[rownames(X_spring) %in% 1967:1988,] 
SST2 <- X_spring[rownames(X_spring) %in% 1989:2012,]
SST3 <- X_spring[rownames(X_spring) %in% 2013:2022,]

CC1 <- CC.FMA[rownames(CC.FMA) %in% 1967:1988, 'stand_bakun_seasonally'] 
CC2 <- CC.FMA[rownames(CC.FMA) %in% 1989:2012, 'stand_bakun_seasonally']
CC3 <- CC.FMA[rownames(CC.FMA) %in% 2013:2022, 'stand_bakun_seasonally']

NCC1 <- NCC.FMA[rownames(NCC.FMA) %in% 1967:1988, 'stand_bakun_seasonally'] 
NCC2 <- NCC.FMA[rownames(NCC.FMA) %in% 1989:2012, 'stand_bakun_seasonally']
NCC3 <- NCC.FMA[rownames(NCC.FMA) %in% 2013:2022, 'stand_bakun_seasonally']

SCC1 <- SCC.FMA[rownames(SCC.FMA) %in% 1967:1988, 'stand_bakun_seasonally'] 
SCC2 <- SCC.FMA[rownames(SCC.FMA) %in% 1989:2012, 'stand_bakun_seasonally']
SCC3 <- SCC.FMA[rownames(SCC.FMA) %in% 2013:2022, 'stand_bakun_seasonally']

GoA1 <- GoA.FMA[rownames(GoA.FMA) %in% 1967:1988, 'stand_bakun_seasonally'] 
GoA2 <- GoA.FMA[rownames(GoA.FMA) %in% 1989:2012, 'stand_bakun_seasonally']
GoA3 <- GoA.FMA[rownames(GoA.FMA) %in% 2013:2022, 'stand_bakun_seasonally']

# calculate separate regressions in each era!

# make objects to catch results
cc.regr1 <- cc.regr2<- cc.regr3 <- ncc.regr1 <- ncc.regr2<- ncc.regr3<-scc.regr1 <- scc.regr2<- scc.regr3<- goa.regr1 <- goa.regr2 <- goa.regr3<- NA

# now loop through each cell

for(i in 1:ncol(SST1)){
  #  i <- 1
  mod <- lm(SST1[,i] ~ CC1)
  cc.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SST2[,i] ~ CC2)
  cc.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SST3[,i] ~ CC3)
  cc.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SST1[,i] ~ NCC1)
  ncc.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST2[,i] ~ NCC2)
  ncc.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST3[,i] ~ NCC3)
  ncc.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST1[,i] ~ SCC1)
  scc.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST2[,i] ~ SCC2)
  scc.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST3[,i] ~ SCC3)
  scc.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST1[,i] ~ GoA1)
  goa.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST2[,i] ~ GoA2)
  goa.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SST3[,i] ~ GoA3)
  goa.regr3[i] <- summary(mod)$coef[2,1]
}
# calculate era differences for each cell

cc.diff <- cc.regr2 - cc.regr1 
ncc.diff <- ncc.regr2 - ncc.regr1 
scc.diff <- scc.regr2 - scc.regr1 
goa.diff <- goa.regr2 - goa.regr1
#diff.lim <- range(pdo.diff, npgo.diff) # limit for plotting

X_cc<- as.data.frame(cc.regr1)%>%
  cbind(as.data.frame(cc.regr2))%>%
  cbind(as.data.frame(cc.regr3))%>%
  cbind(as.data.frame(ncc.regr1))%>%
  cbind(as.data.frame(scc.regr2))%>%
  cbind(as.data.frame(scc.regr3))%>%
  cbind(as.data.frame(scc.regr1))%>%
  cbind(as.data.frame(goa.regr2))%>%
  cbind(as.data.frame(goa.regr3))%>%
  cbind(as.data.frame(goa.regr1))%>%
  cbind(as.data.frame(ncc.regr2))%>%
  cbind(as.data.frame(ncc.regr3))%>%
  #cbind(as.data.frame(cc.diff))%>%
  #cbind(as.data.frame(ncc.diff))%>%
  #cbind(as.data.frame(scc.diff))%>%
  #cbind(as.data.frame(goa.diff))%>%
  rename('CC 1967 - 1988' = cc.regr1, 'NCC 1967 - 1988' = ncc.regr1,'SCC 1967 - 1988' = scc.regr1,  'GoA 1967 - 1988' = goa.regr1,'CC 1989 - 2012' = cc.regr2, 'NCC 1989 - 2012' = ncc.regr2,'SCC 1989 - 2012' = scc.regr2,  'GoA 1989 - 2012' = goa.regr2, 'CC 2013 - 2022' = cc.regr3, 'NCC 2013 - 2022' = ncc.regr3,'SCC 2013 - 2022' = scc.regr3,  'GoA 2013 - 2022' = goa.regr3)

X_cc$latitude <- spring_lat 
X_cc$longitude <- spring_lon+360

X_cc<- X_cc%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")

world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
#plot code
ggplot() + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(120,240), ylim=c(0,60)) +
  geom_raster(data=X_cc, aes(x=longitude,y=latitude,fill = coefficient)) + 
  #facet_wrap(~analysis, ncol = 3) + 
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SST")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 

#### Now SLP ####

# get lat/long
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes
# get times
raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 4:10)] <- NA
dates$winter_year[which(dates$month %in% 11:12)] <- dates$winter_year[which(dates$month %in% 11:12)] + 1 # incremenent 

var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

# The 1st dimension of X is lon, 2nd dimension is lat, 3rd dimension is date. But the array needs to be flattened (time on rows, spatial cells on columns) for easier use.
X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  winter_lat <- lat[-to_drop]
  winter_lon <- lon[-to_drop]
} else {
  winter_lat <- lat
  winter_lon <- lon
}
# Add block for winter average calculations
winter_years <- unique(dates$winter_year)
winter_years <- sort(winter_years[-which(is.na(winter_years))])
X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,] = colMeans(X[which(dates$winter_year == winter_years[i]),])
}
dim(X_winter)

dimnames(X_winter) <- list(as.character(winter_years),paste("N", winter_lat, "E", winter_lon, sep=""))


SLP1 <- X_winter[rownames(X_winter) %in% 1967:1988,] 
SLP2 <- X_winter[rownames(X_winter) %in% 1989:2012,]
SLP3 <- X_winter[rownames(X_winter) %in% 2013:2022,]


# make objects to catch results
cc.regr1 <- cc.regr2<- cc.regr3 <- ncc.regr1 <- ncc.regr2<- ncc.regr3<-scc.regr1 <- scc.regr2<- scc.regr3<- goa.regr1 <- goa.regr2 <- goa.regr3<- NA

# now loop through each cell

for(i in 1:ncol(SLP1)){
  #  i <- 1
  mod <- lm(SLP1[,i] ~ CC1)
  cc.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP2[,i] ~ CC2)
  cc.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP3[,i] ~ CC3)
  cc.regr3[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(SLP1[,i] ~ NCC1)
  ncc.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ NCC2)
  ncc.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ NCC3)
  ncc.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP1[,i] ~ SCC1)
  scc.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ SCC2)
  scc.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ SCC3)
  scc.regr3[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP1[,i] ~ GoA1)
  goa.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP2[,i] ~ GoA2)
  goa.regr2[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(SLP3[,i] ~ GoA3)
  goa.regr3[i] <- summary(mod)$coef[2,1]
}
# calculate era differences for each cell

cc.diff <- cc.regr2 - cc.regr1 
ncc.diff <- ncc.regr2 - ncc.regr1 
scc.diff <- scc.regr2 - scc.regr1 
goa.diff <- goa.regr2 - goa.regr1
#diff.lim <- range(pdo.diff, npgo.diff) # limit for plotting

X_cc_slp<- as.data.frame(cc.regr1)%>%
  cbind(as.data.frame(cc.regr2))%>%
  cbind(as.data.frame(cc.regr3))%>%
  cbind(as.data.frame(ncc.regr1))%>%
  cbind(as.data.frame(scc.regr2))%>%
  cbind(as.data.frame(scc.regr3))%>%
  cbind(as.data.frame(scc.regr1))%>%
  cbind(as.data.frame(goa.regr2))%>%
  cbind(as.data.frame(goa.regr3))%>%
  cbind(as.data.frame(goa.regr1))%>%
  cbind(as.data.frame(ncc.regr2))%>%
  cbind(as.data.frame(ncc.regr3))%>%
  #cbind(as.data.frame(cc.diff))%>%
  #cbind(as.data.frame(ncc.diff))%>%
  #cbind(as.data.frame(scc.diff))%>%
  #cbind(as.data.frame(goa.diff))%>%
  rename('CC 1967 - 1988' = cc.regr1, 'NCC 1967 - 1988' = ncc.regr1,'SCC 1967 - 1988' = scc.regr1,  'GoA 1967 - 1988' = goa.regr1,'CC 1989 - 2012' = cc.regr2, 'NCC 1989 - 2012' = ncc.regr2,'SCC 1989 - 2012' = scc.regr2,  'GoA 1989 - 2012' = goa.regr2, 'CC 2013 - 2022' = cc.regr3, 'NCC 2013 - 2022' = ncc.regr3,'SCC 2013 - 2022' = scc.regr3,  'GoA 2013 - 2022' = goa.regr3)

X_cc_slp$latitude <- winter_lat 
X_cc_slp$longitude <- winter_lon+360

X_cc_slp<- X_cc_slp%>%  
  pivot_longer(!latitude&!longitude,names_to = "analysis", values_to = "coefficient")




world <- st_as_sf(map('world2', plot=F, fill=T))
ggplot() + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(120,240), ylim=c(0,60)) +
  geom_raster(data=X_cc_slp, aes(x=longitude,y=latitude,fill = coefficient)) + 
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  facet_wrap(~analysis, ncol = 3) + 
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("Winter SLP vs. Spring Regional Upwelling")+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(120,240), ylim=c(0,60)) +
  # geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 

##### Data output for analysis summary qmd ####

cor_cc <- X_cc%>%
  mutate(var = 'SST')%>%
  rbind(X_cc_slp%>%mutate(var='SLP'))

saveRDS(cor_cc, file = here('data/physical/correlation_analysis.rds'))
