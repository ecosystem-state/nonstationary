nc <- nc_open(here("som/copernicus_jun27.nc"))

#### Winter SLP ####
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
dim <- c(3,2)

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

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

grid_dim1 <- 3
grid_dim2 <- 2 
grid_shape <- 'rectangular' 

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = scale(X_winter), dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_slp_winter <- data.frame(winter_year = winter_years,
                                          cluster = fit$clustering)%>%
  mutate(var = "SLP_winter")

X_winter <- as.data.frame(X_winter)
X_winter$year <- winter_years

# convert wide to long
winter_long <- pivot_longer(X_winter, cols = 1:(ncol(X_winter)-1))
winter_long$cell_id <- as.numeric(substr(winter_long$name,2,length(winter_long$name)))
winter_long$lat <- winter_lat[winter_long$cell_id]
winter_long$long <- winter_lon[winter_long$cell_id]
winter_long <- dplyr::rename(winter_long, winter_year = year) 
# join in cluster IDs
winter_long <- dplyr::left_join(winter_long, output_slp_winter)

# Now calculated cluster avgs
slp_winter <- dplyr::group_by(winter_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SLP_winter")


#### SLP Winter Globally detrended ####

X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,] = colMeans(X[which(dates$winter_year == winter_years[i]),])
}

# remove annual mean for each year
for(i in 1:nrow(X_winter)) {
  X_winter[i,] = X_winter[i,] - mean(X_winter[i,], na.rm=T)
}
X_winter <- scale(X_winter)

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim


fit <- trainSOM(x.data = X_winter, dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_slp_winter_detrended <- data.frame(winter_year = winter_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SLP_winter_detrended")

X_winter <- as.data.frame(X_winter)
X_winter$year <- winter_years

# convert wide to long
winter_long <- pivot_longer(X_winter, cols = 1:(ncol(X_winter)-1))
winter_long$cell_id <- as.numeric(substr(winter_long$name,2,length(winter_long$name)))
winter_long$lat <- winter_lat[winter_long$cell_id]
winter_long$long <- winter_lon[winter_long$cell_id]
winter_long <- dplyr::rename(winter_long, winter_year = year) 
# join in cluster IDs
winter_long <- dplyr::left_join(winter_long, output_slp_winter_detrended)

# Now calculated cluster avgs
slp_winter_detrended <- dplyr::group_by(winter_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SLP_winter_detrended")


#### Spring SLP ####

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of

raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$spring_year <- dates$year
dates$spring_year[which(dates$month %in% 7:12)] <- NA
dates$spring_year[which(dates$month %in% 1:3)] <- NA

var_name = "msl"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
if(length(to_drop) > 0) {
  X <- X[,-to_drop]# drop cells with NAs 
  spring_lat <- lat[-to_drop]
  spring_lon <- lon[-to_drop]
} else {
  spring_lat <- lat
  spring_lon <- lon
}
# Add block for spring average calculations
spring_years <- unique(dates$spring_year)
spring_years <- sort(spring_years[-which(is.na(spring_years))])
X_spring <- matrix(NA, length(spring_years), ncol(X))
for(i in 1:length(spring_years)) {
  X_spring[i,] = colMeans(X[which(dates$spring_year == spring_years[i]),])
}

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = scale(X_spring), dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_slp_spring <- data.frame(spring_year = spring_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SLP_spring")

X_spring <- as.data.frame(X_spring)
X_spring$year <- spring_years

# convert wide to long
spring_long <- pivot_longer(X_spring, cols = 1:(ncol(X_spring)-1))
spring_long$cell_id <- as.numeric(substr(spring_long$name,2,length(spring_long$name)))
spring_long$lat <- spring_lat[spring_long$cell_id]
spring_long$long <- spring_lon[spring_long$cell_id]
spring_long <- dplyr::rename(spring_long, spring_year = year) 
# join in cluster IDs
spring_long <- dplyr::left_join(spring_long, output_slp_spring)

# Now calculated cluster avgs
slp_spring <- dplyr::group_by(spring_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SLP_spring")


#### SLP Spring Detrended ####

X_spring <- matrix(NA, length(spring_years), ncol(X))
for(i in 1:length(spring_years)) {
  X_spring[i,] = colMeans(X[which(dates$spring_year == spring_years[i]),])
}

# remove annual mean for each year
for(i in 1:nrow(X_spring)) {
  X_spring[i,] = X_spring[i,] - mean(X_spring[i,], na.rm=T)
}
X_spring <- scale(X_spring)

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = X_spring, dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_slp_spring_detrended <- data.frame(spring_year = spring_years,
                                cluster = fit$clustering)%>%
  mutate(var = "SLP_spring_detrended")

X_spring <- as.data.frame(X_spring)
X_spring$year <- spring_years

# convert wide to long
spring_long <- pivot_longer(X_spring, cols = 1:(ncol(X_spring)-1))
spring_long$cell_id <- as.numeric(substr(spring_long$name,2,length(spring_long$name)))
spring_long$lat <- spring_lat[spring_long$cell_id]
spring_long$long <- spring_lon[spring_long$cell_id]
spring_long <- dplyr::rename(spring_long, spring_year = year) 
# join in cluster IDs
spring_long <- dplyr::left_join(spring_long, output_slp_spring_detrended)

# Now calculated cluster avgs
slp_spring_detrended <- dplyr::group_by(spring_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SLP_spring_detrended")


#### Winter SST #### 

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

raw <- ncvar_get(nc, "time")
tunits<-ncatt_get(nc,"time",attname="units")
tustr<-strsplit(tunits$value, " ")
dates <- RNetCDF::utcal.nc("hours since 1900-01-01 00:00:00.0", raw)
dates <- as.data.frame(dates[,c("year","month")])
dates$index <- seq(1, nrow(dates))
dates$winter_year <- dates$year
dates$winter_year[which(dates$month %in% 4:10)] <- NA
dates$winter_year[which(dates$month %in% 11:12)] <- dates$winter_year[which(dates$month %in% 11:12)] + 1 # incremenent 

var_name = "sst"

tmp_array <- ncvar_get(nc,var_name)
dlname <- ncatt_get(nc,var_name,"long_name")
dunits <- ncatt_get(nc,var_name,"units")
fillvalue <- ncatt_get(nc,var_name,"_FillValue")
# replace netCDF fill values with NA's
tmp_array[tmp_array==fillvalue$value] <- NA
X <- tmp_array

X <- aperm(X, 3:1) # transpose array # X
X <- matrix(X, nrow=dim(X)[1], ncol=prod(dim(X)[2:3])) # months in rows, cells in columns! 
to_drop <- which(is.na(apply(X,2,sum)))
X <- X[,-to_drop]# drop cells with NAs 
winter_lat <- lat[-to_drop]
winter_lon <- lon[-to_drop]
# Add block for winter average calculations
winter_years <- unique(dates$winter_year)
winter_years <- sort(winter_years[-which(is.na(winter_years))])
X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,] = colMeans(X[which(dates$winter_year == winter_years[i]),])
}

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = scale(X_winter), dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_sst_winter <- data.frame(winter_year = winter_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SST_winter")

X_winter <- as.data.frame(X_winter)
X_winter$year <- winter_years

# convert wide to long
winter_long <- pivot_longer(X_winter, cols = 1:(ncol(X_winter)-1))
winter_long$cell_id <- as.numeric(substr(winter_long$name,2,length(winter_long$name)))
winter_long$lat <- winter_lat[winter_long$cell_id]
winter_long$long <- winter_lon[winter_long$cell_id]
winter_long <- dplyr::rename(winter_long, winter_year = year) 
# join in cluster IDs
winter_long <- dplyr::left_join(winter_long, output_sst_winter)

# Now calculated cluster avgs
sst_winter <- dplyr::group_by(winter_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SST_winter")

#### SST Winter Detrended ####

X_winter <- matrix(NA, length(winter_years), ncol(X))
for(i in 1:length(winter_years)) {
  X_winter[i,] = colMeans(X[which(dates$winter_year == winter_years[i]),])
}

# remove annual mean for each year
for(i in 1:nrow(X_winter)) {
  X_winter[i,] = X_winter[i,] - mean(X_winter[i,], na.rm=T)
}
X_winter <- scale(X_winter)

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = X_winter, dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_sst_winter_detrended <- data.frame(winter_year = winter_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SST_winter_detrended")

X_winter <- as.data.frame(X_winter)
X_winter$year <- winter_years

# convert wide to long
winter_long <- pivot_longer(X_winter, cols = 1:(ncol(X_winter)-1))
winter_long$cell_id <- as.numeric(substr(winter_long$name,2,length(winter_long$name)))
winter_long$lat <- winter_lat[winter_long$cell_id]
winter_long$long <- winter_lon[winter_long$cell_id]
winter_long <- dplyr::rename(winter_long, winter_year = year) 
# join in cluster IDs
winter_long <- dplyr::left_join(winter_long, output_sst_winter_detrended)

# Now calculated cluster avgs
sst_winter_detrended <-  dplyr::group_by(winter_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SST_winter_detrended")

#### Spring SST ####

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

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

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = scale(X_spring), dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_sst_spring <- data.frame(spring_year = spring_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SST_spring")

X_spring <- as.data.frame(X_spring)
X_spring$year <- spring_years

# convert wide to long
spring_long <- pivot_longer(X_spring, cols = 1:(ncol(X_spring)-1))
spring_long$cell_id <- as.numeric(substr(spring_long$name,2,length(spring_long$name)))
spring_long$lat <- spring_lat[spring_long$cell_id]
spring_long$long <- spring_lon[spring_long$cell_id]
spring_long <- dplyr::rename(spring_long, spring_year = year) 
# join in cluster IDs
spring_long <- dplyr::left_join(spring_long, output_sst_spring)

# Now calculated cluster avgs
sst_spring<- dplyr::group_by(spring_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SST_spring")


#### SST Spring Detrended ####

X_spring <- matrix(NA, length(spring_years), ncol(X))
for(i in 1:length(spring_years)) {
  X_spring[i,] = colMeans(X[which(dates$spring_year == spring_years[i]),])
}

# remove annual mean for each year
for(i in 1:nrow(X_spring)) {
  X_spring[i,] = X_spring[i,] - mean(X_spring[i,], na.rm=T)
}
X_spring <- scale(X_spring)

grids <- expand.grid(xdim = 2:7)
grids$ydim <- grids$xdim + 1
grids <- rbind(grids, data.frame(xdim = 2:7, ydim = 2:7))
grids$mean_dist <- NA
grids$size <- grids$xdim * grids$ydim

fit <- trainSOM(x.data = X_spring, dimension = dim, verbose = FALSE, nb.save = 0, topo = "square")

output_sst_spring_detrended <- data.frame(spring_year = spring_years,
                     cluster = fit$clustering)%>%
  mutate(var = "SST_spring_detrended")

X_spring <- as.data.frame(X_spring)
X_spring$year <- spring_years

# convert wide to long
spring_long <- pivot_longer(X_spring, cols = 1:(ncol(X_spring)-1))
spring_long$cell_id <- as.numeric(substr(spring_long$name,2,length(spring_long$name)))
spring_long$lat <- spring_lat[spring_long$cell_id]
spring_long$long <- spring_lon[spring_long$cell_id]
spring_long <- dplyr::rename(spring_long, spring_year = year) 
# join in cluster IDs
spring_long <- dplyr::left_join(spring_long, output_sst_spring_detrended)

# Now calculated cluster avgs
sst_spring_detrended  <- dplyr::group_by(spring_long, cell_id) %>%
  dplyr::mutate(z = (value - mean(value))/sd(value)) %>%
  dplyr::group_by(cell_id, cluster) %>%
  dplyr::summarise(mean = mean(z),
                   lat = lat[1],
                   long = long[1])%>%
  mutate(var = "SST_spring_detrended")

##### Data output for analysis summary qmd ####

SOM_output <- output_slp_winter%>%rename(year=winter_year)%>%
  rbind(output_slp_winter_detrended%>%rename(year=winter_year))%>%
  rbind(output_slp_spring%>%rename(year=spring_year))%>%
  rbind(output_slp_spring_detrended%>%rename(year=spring_year))%>%
  rbind(output_sst_winter%>%rename(year=winter_year))%>%
  rbind(output_sst_winter_detrended%>%rename(year=winter_year))%>%
  rbind(output_sst_spring%>%rename(year=spring_year))%>%
  rbind(output_sst_spring_detrended%>%rename(year=spring_year))

saveRDS(SOM_output, file = here('data/physical/SOM_output.rds'))

SOM_data <- slp_winter%>%
  rbind(slp_winter_detrended)%>%
  rbind(slp_spring)%>%
  rbind(slp_spring_detrended)%>%
  rbind(sst_winter)%>%
  rbind(sst_winter_detrended)%>%
  rbind(sst_spring)%>%
  rbind(sst_spring_detrended)

saveRDS(SOM_data, file = here('data/physical/SOM_data.rds'))

