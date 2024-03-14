library(terra) 

ddd <-readRDS('data/physical/correlation_analysis_diff.rds')

world <- st_as_sf(map('world2', plot=F, fill=T)) #base layer for land masses
world3<-world%>%filter(ID=="Canada"|ID=="USA"|ID=="Mexico") #pulling ou USA,M and C
world3<-st_union(world3) #making them one multipoly object
buff<-buffer(vect(world3),width=130000) #adding a buffer
buff2<- st_as_sf(buff) #making it a sf object

buffer<- data.frame(st_coordinates(buff2))%>% #creating the buffer polygon
  filter(Y>30&Y<50)%>%
  mutate(X=X+360)%>%
  filter((X>220&X<250))%>%
 # filter(L1==1,L2==1,L3==1020)%>%
 # select(X,Y)%>%
  mutate(x=X,y=Y)
unique(coast$L1)
unique(coast$L2)
unique(coast$L3)
coast2<-st_coordinates(st_as_sf(buffer(vect(world3),width=1)))
coast<- data.frame(coast2)%>%
  filter(Y>30&Y<55)%>%
  mutate(X=X+360)%>%
  filter((X>220&X<250))%>%
  filter((L2==113))%>%
 # filter(L1==1,L2==1,L3==1020)%>%|L2==151 OR56
 # select(X,Y)%>%
  purrr::map_df(rev)%>%
  mutate(x=X,y=Y)

#no 31
buffer<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>31&y<48.5)
bufferNCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>40.4401&y<48.5)
bufferCCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>34.4486&y<40.5401)
bufferSCC<-rbind(buffer[8:228,],coast[1:2000,],buffer[8,])%>%
  filter(y>31&y<34.5486)
analysis2<-data.frame(analysis=unique(X_PDO_SLP$analysis), analysis2=c("1967 - 1988","1989 - 2012","2013 - 2023","X",
                                                            "1967 - 1988","1989 - 2012","2013 - 2023", "X",
                                                            "1967 - 1988", "1989 - 2012","2013 - 2023","X",
                                                            "1967 - 1988","1989 - 2012","2013 - 2023","X"),
                      Index=c(rep('PDO', 4), rep('ONI',4), rep('NPGO',4), rep('NPH',4))) 
analysis2<-merge(X_PDO_SLP,analysis2)%>%
  filter(analysis2!='X')
SLP_coast2<- ggplot() + 
 geom_raster(data=analysis2, aes(x=longitude,y=latitude,fill = coefficient)) + 
   facet_grid(Index~analysis2) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(230,245), ylim=c(30,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("Spring SLP Anomalies (Pa) vs. Climate Indices")+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks = c(230,240))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(strip.background=element_rect(colour="black",
                                    fill="white"),panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))

pdf(file = "Output/Figures/SLP_coast2.pdf",   # The directory you want to save the file in
    width = 7, # The width of the plot in inches
    height =7)
SLP_coast2
dev.off()

SLP_coast<- ggplot() + 
 geom_raster(data=X_PDO_SLP, aes(x=longitude,y=latitude,fill = coefficient)) + 
   facet_wrap(~analysis, ncol = 4) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(233,245), ylim=c(31,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
    ggtitle("Spring SLP Anomalies (Pa) vs. Climate Indices")+
  ylab("")+
  xlab("")+
  scale_x_continuous(breaks = c(125, 126))+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA))
SLP_coast



SST_coast<- ggplot() + 
 geom_raster(data=X_PDO_SST, aes(x=longitude,y=latitude,fill = coefficient)) + 
   facet_wrap(~analysis, ncol = 3) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(233,245), ylim=c(31,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SST")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SST_coast

SLP_diff_coast<- ggplot() + 
 geom_raster(data=X_Diff_SLP, aes(x=longitude,y=latitude,fill = coefficient)) + 
   facet_wrap(~analysis, ncol = 3) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_sf(data=world, col="black", fill="darkgoldenrod3") +
  coord_sf(xlim=c(233,245), ylim=c(31,50)) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  ggtitle("SLP")+
  #geom_contour(data=X_cc, aes(x=longitude,y=latitude,z = coefficient), col="lightgrey", lwd=0.5)+
  theme(plot.title = element_text(hjust = 0.5), panel.border = element_rect(fill = NA)) 
SLP_diff_coast

# Coastal ocean subregion 1
regional.coefs<-NA
unique(X_Diff_SLP$analysis)
ncc_sf <- sf_polygon(bufferNCC%>%dplyr::select(X,Y))
ccc_sf <- sf_polygon(bufferCCC%>%dplyr::select(X,Y))
scc_sf <- sf_polygon(bufferSCC%>%dplyr::select(X,Y))

regional_coef_diff<-function(regional_sf) {
for(i in 1:length(unique(X_Diff_SLP$analysis))){
  time<-data.frame(unique(X_Diff_SLP$analysis))[i,1]
  X<-data.frame(X_Diff_SLP%>%filter(analysis==time)%>%
     dplyr::select(longitude,latitude,  coefficient))
  r<-st_as_sf(X,coords = c("longitude","latitude"))
  coef_in_tract <- st_join(r,regional_sf,join = st_within)
  df<-na.omit(coef_in_tract)
  plot(df)
  df.data<-data.frame(cbind(coefficient=data.frame(df)$coefficient,st_coordinates(df),
                          analysis=time))
colnames(df.data)<-c("coef",'X',"Y") 
regional.coefs<-rbind(regional.coefs, df.data)
}
  return(regional.coefs)
}

ncc_coef_diff<-na.omit(cbind(regional_coef_diff(ncc_sf),region="NCC"))
colnames(ncc_coef_diff)<-c("coef","X","Y","analysis","region")
ccc_coef_diff<-na.omit(cbind(regional_coef_diff(ccc_sf),region="CCC"))
colnames(ccc_coef_diff)<-c("coef","X","Y","analysis","region")
scc_coef_diff<-na.omit(cbind(regional_coef_diff(scc_sf),region="SCC"))
colnames(scc_coef_diff)<-c("coef","X","Y","analysis","region")

coef_diff<-ncc_coef_diff%>%add_row(ccc_coef_diff)%>%add_row(scc_coef_diff)

diff_summary<- coef_diff%>%dplyr::select(-X, -Y)%>%group_by( region,analysis)%>%
  summarise(mean=mean(as.numeric(coef)), sd=sd(as.numeric(coef)))

SLP_diff_coast<-ggplot() + 
  geom_raster(data=coef_diff, aes(x=X,y=Y,fill = coef)) + 
  facet_wrap(~analysis, ncol = 4) +
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)

ncc_coef_diff



ncc_sf <- sf_polygon(ncc)
ncc<-bufferNCC%>%rename(longitude=X, latitude=Y)%>%dplyr::select(longitude,latitude)

X<-data.frame(X_Diff_SLP%>%filter(analysis=="F. ONI 2013:2023 - 1967:2023")%>%
         dplyr::select(longitude,latitude,  coefficient))

r<-st_as_sf(X,coords = c("longitude","latitude"))

tree_in_tract <- st_join(r,ncc_sf,join = st_within)
df<-na.omit(tree_in_tract)
plot(tree_in_tract)
df.data<-data.frame(cbind(coefficient=data.frame(df)$coefficient,st_coordinates(df)))
colnames(df.data)<-c("coef",'X',"Y")
ggplot() + 
  geom_raster(data=df.data, aes(x=X,y=Y,fill = coef)) + 
  geom_polygon(data=bufferNCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferCCC,aes(x=x,y=y),color='black',fill=NA)+
  geom_polygon(data=bufferSCC,aes(x=x,y=y),color='black',fill=NA)
