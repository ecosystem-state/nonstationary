


climate_dat_cop<-readRDS(here('data/physical/climate_dat_cop.rds'))
dfa<-readRDS(here('data/physical/climate_dat_dfa.rds'))%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))
dat <- dfa%>%dplyr::select(estimate, Year_lag, lower, upper,trend,season)%>%
  mutate(Year_lag = Year_lag+1)%>%
  rename(estimateoffset1=estimate, loweroffset1=lower, upperoffset1=upper)%>%
  left_join(dfa, by=c('trend', 'Year_lag','season'))

data_spring <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
         seasonal_copepod_southern)%>%
  left_join(climate_dat_cop%>%filter(season=="Spring")%>%
              dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
              distinct()%>%
              filter(Year_lag<2023))%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern,
         seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  pivot_longer(!c(Year_lag,  NPH, NPGO,PDO,ONI), 
               names_to = "trend", values_to = "estimate")%>%
  mutate(season="Spring")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

data_winter <- climate_dat_cop%>%filter(season=="Summer")%>%
  dplyr::select(Year_lag, region, seasonal_copepod_northern, 
                seasonal_copepod_southern)%>%
  left_join(climate_dat_cop%>%filter(season=="Winter")%>%
              dplyr::select(Year_lag, region, seasonal_PDO,  seasonal_NPGO,  seasonal_ONI,  seasonal_NPH)%>%
              distinct()%>%
              filter(Year_lag<2023))%>%
  distinct()%>%
  mutate(period=ifelse(Year_lag<2013,2,3))%>%
  filter(Year_lag>1996)%>%
  dplyr::select(Year_lag, seasonal_copepod_northern,seasonal_copepod_southern,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  pivot_longer(!c(Year_lag,  NPH, NPGO,PDO,ONI), 
               names_to = "trend", values_to = "estimate")%>%
  mutate(season="Winter")%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value")


dat.long.offset<- dat%>%
  dplyr::select(Year_lag, season, estimateoffset1,trend,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimateoffset1,trend), 
               names_to = "Index_Name", values_to = "Index_Value")

##### Spring Biological Data Plot #####

dat.long<- bind_rows(dat%>%filter(season=="Spring")%>%
  dplyr::select(Year_lag, season, estimate,trend,
                seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
  rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
  distinct()%>%
  pivot_longer(!c(Year_lag, season, estimate,trend), 
               names_to = "Index_Name", values_to = "Index_Value"),
  dat%>%filter(season=="Winter")%>%
    dplyr::select(Year_lag, season, estimate,trend,
                  seasonal_NPH,seasonal_NPGO,seasonal_PDO,seasonal_ONI)%>%
    rename(NPH=seasonal_NPH,NPGO=seasonal_NPGO,PDO=seasonal_PDO,ONI=seasonal_ONI)%>%
    distinct()%>%
    pivot_longer(!c(Year_lag, season, estimate,trend), 
                 names_to = "Index_Name", values_to = "Index_Value"))

dat.lm<-data_spring%>%
  bind_rows(dat.long,data_winter)%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  filter(trend!="SEA")%>%
  mutate(trend=ifelse(trend=="CALCOFI","CALCOFI (SCC)",
                      ifelse(trend=="RREAS", "RREAS (CCC)", 
                             ifelse(trend=="seasonal_copepod_southern","Southern Copepod (NCC)","Northern Copepod (NCC)"))))%>%
  mutate(trend = fct_relevel(trend, "Northern Copepod (NCC)", 
                             "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))


survey.lm<-ggplot(data = dat.lm%>%filter(season=='Spring'), aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
survey.lm


##### Winter Biological Data Plot #####


survey.lm.w<-ggplot(data = dat.lm%>%filter(season=='Winter'), aes(y = estimate, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
survey.lm.w

##### Spring Biological Data Plot #####

dat.long.offset<- dat.long.offset%>%
  mutate(period=ifelse(Year_lag<1989,1,ifelse(Year_lag>2012,3,2)))%>%
  filter(trend!="SEA")%>%
  mutate(trend=ifelse(trend=="CALCOFI","CALCOFI (SCC)",
                      ifelse(trend=="RREAS", "RREAS (CCC)", 
                             ifelse(trend=="seasonal_copepod_southern",
                                    "Southern Copepod (NCC)","Northern Copepod (NCC)"))))%>%
  mutate(trend = fct_relevel(trend, "Northern Copepod (NCC)", 
                             "Southern Copepod (NCC)","RREAS (CCC)","CALCOFI (SCC)"))

survey.lm.offset<-ggplot(data = dat.lm.off%>%filter(season=='Spring'), aes(y = estimateoffset1, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~trend, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Index of Abundance") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-year lag")
survey.lm.offset


##### TUMI  ####
climate_dat <-readRDS(here('data/physical/climate_dat_upwelling.rds'))
data_up<- climate_dat%>%filter(region!="GoA")%>%
  select(!c(stand_bakun_seasonally,era.region,annual_NPGO, annual_NPH, annual_ONI))%>%
  rename(NPGO=seasonal_NPGO, ONI=seasonal_ONI, NPH=seasonal_NPH, PDO=seasonal_PDO)%>%
  pivot_longer(!c(Year_lag, region, season, period, stand_tumi, stand_sti,stand_lusi), 
               names_to = "Index_Name", values_to = "Index_Value")%>%
  rename(TUMI=stand_tumi, STI=stand_sti, LUSI=stand_lusi)%>%
  pivot_longer(!c(Year_lag, region, season, period, Index_Value,Index_Name), 
               names_to = "Upwelling", values_to = "Upwelling_Index")%>%
  mutate(region = fct_relevel(region, "Northern CC", "Central CC", "Southern CC"))%>%
  distinct()

tumi.spring<-ggplot(data = data_up%>%filter(season=='Spring'&Upwelling=="TUMI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "TUMI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
tumi.spring

tumi.winter<-ggplot(data = data_up%>%filter(season=='Winter'&Upwelling=="TUMI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "TUMI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
tumi.winter

##### STI  ####

sti.spring<-ggplot(data = data_up%>%filter(season=='Spring'&Upwelling=="STI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "STI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("STI Spring")
sti.spring

sti.winter<-ggplot(data = data_up%>%filter(season=='Winter'&Upwelling=="STI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "STI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Winter")
sti.winter


##### LUSI  ####

lusi.spring<-ggplot(data = data_up%>%filter(season=='Spring'&Upwelling=="LUSI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "LUSI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
lusi.spring

lusi.winter<-ggplot(data = data_up%>%filter(season=='Winter'&Upwelling=="LUSI"), aes(y = Upwelling_Index, x =Index_Value,col=as.factor(period))) +
  facet_grid(Index_Name~region, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "LUSI") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Climate Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("LUSI")
lusi.winter

##### Biological Upwelling #####

bioup<-left_join(dat.lm%>%
                   select(!c(period,Index_Name, Index_Value))%>%distinct()%>%
                   mutate(region=ifelse(trend=='CALCOFI','Southern CC',
                                        ifelse(trend=="RREAS", 'Central CC','Northern CC'))), 
                 data_up%>%select(!c(Index_Name,Index_Value))%>%
                   distinct()%>%
                   pivot_wider(names_from = Upwelling, values_from = Upwelling_Index))%>%
  pivot_longer(!c(Year_lag, trend,estimate,region, season, period), 
               names_to = "Upwelling", values_to = "Upwelling_Index")

bioup_plot<-ggplot(data = bioup%>%
                     filter(season=='Spring'), aes(y = estimate, x =Upwelling_Index,col=as.factor(period))) +
  facet_grid(trend~Upwelling, scales='free') +
  geom_point(aes(col=as.factor(period))) +
  # geom_text(aes(label=Year_lag,col=as.factor(period))) +
  geom_smooth(method = "lm", se = FALSE, aes(col=as.factor(period))) +
  # geom_smooth(method = "lm", se = FALSE, col='grey') +
  scale_y_continuous(name = "Biological Index Value") +
  scale_color_manual(values =  col[1:3], name="Period",labels=c('1967 - 1988', '1989 - 2012','2013 - 2022'))+
  theme_bw()+
  xlab("Upwelling Index Value")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("Spring")
bioup_plot


##### Beta Plots #####

Violin_Data <-readRDS('data/Violin_Data.rds')%>%
  mutate(region = fct_relevel(region, "Northern CC", "Central CC", "Southern CC"))%>%
  mutate(survey = fct_relevel(survey, "N. Copepod", 
                             "S. Copepod","RREAS","CALCOFI", 
                           "TUMI","STI","LUSI","Upwelling"))

ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="TUMI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior Density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("TUMI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="TUMI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("TUMI Winter")

ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="STI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  ggtitle("STI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="STI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  ggtitle("STI Winter")


ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="LUSI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("LUSI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="LUSI"), aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("LUSI Winter")


ggplot(Violin_Data%>%filter(Season=="Spring"&period!=4)%>%
                    filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
                    filter(Index=="PDO"|Index=="NPH"|Index=="NPGO"|Index=="ONI"), 
       aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")

ggplot(Violin_Data%>%filter(Season=="Winter"&period!=4)%>%
         filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="PDO"|Index=="NPH"|Index=="NPGO"|Index=="ONI"), 
       aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")

ggplot(Violin_Data%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="TUMI"|Index=="LUSI"|Index=="STI"), 
       aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")



ggplot(Violin_Data%>%filter(Season=="Spring"&lag==1&period!=4), 
       aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (beta)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-year lag")

##### Alpha Plots #####
Violin_Data <-readRDS('data/Violin_Data.rds')%>%
  mutate(region = fct_relevel(region, "Northern CC", "Central CC", "Southern CC"))%>%
  mutate(survey = fct_relevel(survey, "N. Copepod", 
                              "S. Copepod","RREAS","CALCOFI", 
                              "TUMI","STI","LUSI","Upwelling"))

ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="TUMI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior Density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("TUMI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="TUMI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("TUMI Winter")

ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="STI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Intercept (alpha)",
       y = "Posterior density")+
  ggtitle("STI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="STI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  theme(plot.title = element_text(hjust = 0.5))+
  labs(x = "Intercept (alpha)",
       y = "Posterior density")+
  ggtitle("STI Winter")


ggplot(Violin_Data%>%filter(Season=="Spring"&survey=="LUSI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("LUSI Spring")

ggplot(Violin_Data%>%filter(Season=="Winter"&survey=="LUSI"), aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~region) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Inetercept (alpha)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("LUSI Winter")


ggplot(Violin_Data%>%filter(Season=="Spring"&period!=4)%>%
         filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="PDO"|Index=="NPH"|Index=="NPGO"|Index=="ONI"), 
       aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")

ggplot(Violin_Data%>%filter(Season=="Winter"&period!=4)%>%
         filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="PDO"|Index=="NPH"|Index=="NPGO"|Index=="ONI"), 
       aes(x = beta,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")

ggplot(Violin_Data%>%filter(survey=="CALCOFI"|survey=="RREAS"|survey=="N. Copepod"|survey=="S. Copepod")%>%
         filter(Index=="TUMI"|Index=="LUSI"|Index=="STI"), 
       aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")

ggplot(Violin_Data%>%filter(Season=="Spring"&lag==1&period!=4), 
       aes(x = alpha,fill = as.factor(period), group=as.factor(period))) +
  theme_bw() +
  facet_grid(Index~survey, scales='free') +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(col[1],col[2],col[3],'grey')) +
  #theme(legend.title = element_blank(), legend.position = 'top', legend.key.size = unit(3, 'mm')) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (alpha)",
       y = "Posterior density")+
  theme(plot.title = element_text(hjust = 0.5))+
  ggtitle("1-year lag")


