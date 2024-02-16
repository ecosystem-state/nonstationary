library(here)
library(dplyr)
library(tidyverse)
library(rerddapXtracto)


# Coastal ocean subregion 1: Puget Sound 
poly1 <- read_csv(here("polygons_lat_lon", "poly1_coord.csv"))

# Coastal ocean subregion 2: WA Coast & Columbia R
poly2 <- read_csv(here("polygons_lat_lon", "poly2_coord.csv"))

# Coastal ocean subregion 3: N and Central OR Coast
poly3 <- read_csv(here("polygons_lat_lon", "poly3_coord.csv"))

# Coastal ocean subregion 4: S OR and N CA Coast
poly4 <- read_csv(here("polygons_lat_lon", "poly4_coord.csv"))

# Coastal ocean subregion 5: Mendocino Coast and San Fran Coast
poly5 <- read_csv(here("polygons_lat_lon", "poly5_coord.csv"))




# ERSST

# NOAA ERSSTv5 (in situ only), 2°, Global, Monthly, 1854-present 
dataInfo <- rerddap::info('nceiErsstv5')
tpos<-c('1854-01-01','2023-11-15')  


# Coastal ocean subregion 1
ypos<-(poly1[,2])       # lat
xpos<-(poly1[,1]) +360  # lon + 360
zpos<-0                 # depth

ERSST_poly1 <- rxtracto_3D(dataInfo, parameter='sst', xcoord=xpos, ycoord=ypos, zcoord=zpos, tcoord=tpos, # sst data in array of [Lon, Lat, Time]
                           xName='longitude', yName='latitude', zName='depth', tName='time', 
                           verbose=TRUE) 
ERSST_poly1_yr_mo_df <- data.frame(year = as.numeric(format(ERSST_poly1$time, format="%Y")), 
                           month = as.numeric(format(ERSST_poly1$time, format="%m")),
                           sst = apply(ERSST_poly1$sst, 4, mean, na.rm = TRUE) )

ERSST_poly1_yr_mo_df <-ERSST_poly1_yr_mo_df%>%add_column(ecoregion=1)


# Coastal ocean subregion 2
ypos<-(poly2[,2])       # lat
xpos<-(poly2[,1]) +360  # lon + 360
zpos<-0                 # depth

ERSST_poly2 <- rxtracto_3D(dataInfo, parameter='sst', xcoord=xpos, ycoord=ypos, zcoord=zpos, tcoord=tpos, # sst data in array of [Lon, Lat, Time]
                           xName='longitude', yName='latitude', zName='depth', tName='time', 
                           verbose=TRUE) 
ERSST_poly2_yr_mo_df <- data.frame(year = as.numeric(format(ERSST_poly2$time, format="%Y")), 
                           month = as.numeric(format(ERSST_poly2$time, format="%m")),
                           sst = apply(ERSST_poly2$sst, 4, mean, na.rm = TRUE) )

ERSST_poly2_yr_mo_df <-ERSST_poly2_yr_mo_df%>%add_column(ecoregion=2)

# Coastal ocean subregion 3
ypos<-(poly3[,2])       # lat
xpos<-(poly3[,1]) +360  # lon + 360
zpos<-0                 # depth

ERSST_poly3 <- rxtracto_3D(dataInfo, parameter='sst', xcoord=xpos, ycoord=ypos, zcoord=zpos, tcoord=tpos, # sst data in array of [Lon, Lat, Time]
                           xName='longitude', yName='latitude', zName='depth', tName='time', 
                           verbose=TRUE) 
ERSST_poly3_yr_mo_df <- data.frame(year = as.numeric(format(ERSST_poly3$time, format="%Y")), 
                                   month = as.numeric(format(ERSST_poly3$time, format="%m")),
                                   sst = apply(ERSST_poly3$sst, 4, mean, na.rm = TRUE) )

ERSST_poly3_yr_mo_df <-ERSST_poly3_yr_mo_df%>%add_column(ecoregion=3)



# Coastal ocean subregion 4
ypos<-(poly4[,2])       # lat
xpos<-(poly4[,1]) +360  # lon + 360
zpos<-0                 # depth

ERSST_poly4 <- rxtracto_3D(dataInfo, parameter='sst', xcoord=xpos, ycoord=ypos, zcoord=zpos, tcoord=tpos, # sst data in array of [Lon, Lat, Time]
                           xName='longitude', yName='latitude', zName='depth', tName='time', 
                           verbose=TRUE) 
ERSST_poly4_yr_mo_df <- data.frame(year = as.numeric(format(ERSST_poly4$time, format="%Y")), 
                                   month = as.numeric(format(ERSST_poly4$time, format="%m")),
                                   sst = apply(ERSST_poly4$sst, 4, mean, na.rm = TRUE) )

ERSST_poly4_yr_mo_df <-ERSST_poly4_yr_mo_df%>%add_column(ecoregion=4)

# Coastal ocean subregion 5
ypos<-(poly5[,2])       # lat
xpos<-(poly5[,1]) +360  # lon + 360
zpos<-0                 # depth

ERSST_poly5 <- rxtracto_3D(dataInfo, parameter='sst', xcoord=xpos, ycoord=ypos, zcoord=zpos, tcoord=tpos, # sst data in array of [Lon, Lat, Time]
                           xName='longitude', yName='latitude', zName='depth', tName='time', 
                           verbose=TRUE) 
ERSST_poly5_yr_mo_df <- data.frame(year = as.numeric(format(ERSST_poly5$time, format="%Y")), 
                                   month = as.numeric(format(ERSST_poly5$time, format="%m")),
                                   sst = apply(ERSST_poly5$sst, 4, mean, na.rm = TRUE) )
ERSST_poly5_yr_mo_df <-ERSST_poly5_yr_mo_df%>%add_column(ecoregion=5)

ERSST_poly<-ERSST_poly1_yr_mo_df%>%add_row(ERSST_poly2_yr_mo_df)%>%
                                           add_row(ERSST_poly3_yr_mo_df)%>%
                                           add_row(ERSST_poly4_yr_mo_df)%>%
                                           add_row(ERSST_poly5_yr_mo_df)

saveRDS(ERSST_poly, file = "data/environment/SST/ERSST_poly.rds")

# ----------- satellite SST -----------

# Multi-scale Ultra-high Resolution (MUR) SST Analysis fv04.1, Global, 0.01°, 2002-present, Monthly
dataInfo <- rerddap::info('jplMURSST41mday')
tpos<-c('2003-01-01','2023-12-16')  


# Coastal ocean subregion 1
ypos<-(poly1[,2])    # lat
xpos<-(poly1[,1])    # lon

SST_MUR_poly1 <- rxtractogon(dataInfo, xcoord=xpos, ycoord=ypos, tcoord=tpos, parameter='sst') # sst data in array of [Lon, Lat, Time]
SST_MUR_poly1_yr_mo_df <- data.frame(year = as.numeric(format(SST_MUR_poly1$time, format="%Y")), 
                                     month = as.numeric(format(SST_MUR_poly1$time, format="%m")),
                                     sst = apply(SST_MUR_poly1$sst, 3, mean, na.rm = TRUE) )
SST_MUR_poly1_yr_mo_df <-SST_MUR_poly1_yr_mo_df%>%add_column(ecoregion=5)

# Coastal ocean subregion 2
ypos<-(poly2[,2])    # lat
xpos<-(poly2[,1])    # lon  

SST_MUR_poly2 <- rxtractogon(dataInfo, xcoord=xpos, ycoord=ypos, tcoord=tpos, parameter='sst') # sst data in array of [Lon, Lat, Time]
SST_MUR_poly2_yr_mo_df <- data.frame(year = as.numeric(format(SST_MUR_poly2$time, format="%Y")), 
                                     month = as.numeric(format(SST_MUR_poly2$time, format="%m")),
                                     sst = apply(SST_MUR_poly2$sst, 3, mean, na.rm = TRUE) )
SST_MUR_poly2_yr_mo_df <-SST_MUR_poly2_yr_mo_df%>%add_column(ecoregion=5)


# Coastal ocean subregion 3
ypos<-(poly3[,2])    # lat
xpos<-(poly3[,1])    # lon

SST_MUR_poly3 <- rxtractogon(dataInfo, xcoord=xpos, ycoord=ypos, tcoord=tpos, parameter='sst') # sst data in array of [Lon, Lat, Time]
SST_MUR_poly3_yr_mo_df <- data.frame(year = as.numeric(format(SST_MUR_poly3$time, format="%Y")), 
                                     month = as.numeric(format(SST_MUR_poly3$time, format="%m")),
                                     sst = apply(SST_MUR_poly3$sst, 3, mean, na.rm = TRUE) )

SST_MUR_poly3_yr_mo_df <-SST_MUR_poly3_yr_mo_df%>%add_column(ecoregion=5)

# Coastal ocean subregion 4
ypos<-(poly4[,2])    # lat
xpos<-(poly4[,1])    # lon

SST_MUR_poly4 <- rxtractogon(dataInfo, xcoord=xpos, ycoord=ypos, tcoord=tpos, parameter='sst') # sst data in array of [Lon, Lat, Time]
SST_MUR_poly4_yr_mo_df <- data.frame(year = as.numeric(format(SST_MUR_poly4$time, format="%Y")), 
                                     month = as.numeric(format(SST_MUR_poly4$time, format="%m")),
                                     sst = apply(SST_MUR_poly4$sst, 3, mean, na.rm = TRUE) )

SST_MUR_poly4_yr_mo_df <-SST_MUR_poly4_yr_mo_df%>%add_column(ecoregion=5)

# Coastal ocean subregion 5
ypos<-(poly5[,2])    # lat
xpos<-(poly5[,1])    # lon 

SST_MUR_poly5 <- rxtractogon(dataInfo, xcoord=xpos, ycoord=ypos, tcoord=tpos, parameter='sst') # sst data in array of [Lon, Lat, Time]
SST_MUR_poly5_yr_mo_df <- data.frame(year = as.numeric(format(SST_MUR_poly5$time, format="%Y")), 
                                     month = as.numeric(format(SST_MUR_poly5$time, format="%m")),
                                     sst = apply(SST_MUR_poly5$sst, 3, mean, na.rm = TRUE) )


SST_MUR_poly5_yr_mo_df <-SST_MUR_poly5_yr_mo_df%>%add_column(ecoregion=5)

SST_MUR_poly<-SST_MUR_poly1_yr_mo_df%>%add_row(SST_MUR_poly2_yr_mo_df)%>%
                                           add_row(SST_MUR_poly3_yr_mo_df)%>%
                                           add_row(SST_MUR_poly4_yr_mo_df)%>%
                                           add_row(SST_MUR_poly5_yr_mo_df)

saveRDS(SST_MUR_poly, file = "data/environment/SST/SST_MUR_poly.rds")

