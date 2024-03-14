
Violin_Data <-readRDS(here('data/Violin_Data.rds'))
unique(Violin_Data$survey)
Violin_indices <- filter(Violin_Data, survey!="Upwelling") 
Violin_upwelling <- filter(Violin_Data, survey=="Upwelling") 


ratio_biological<-rbind(Violin_indices%>%filter(period==2)%>%group_by(Index,survey, lag, Season)%>%
                       add_column((Violin_indices%>%filter(period==3)%>%
                                     group_by(Index,survey, lag, Season))$beta-(Violin_indices%>%
                                            filter(period==2)%>%
                                            group_by(Index,survey, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%
  mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2")%>%
  bind_rows(rbind(Violin_indices%>%filter(period==3)%>%group_by(Index,region, lag, Season)%>%
                       add_column((Violin_indices%>%filter(period==3)%>%
                                     group_by(Index,region, lag, Season))$beta-(Violin_indices%>%
                                            filter(period==4)%>%
                                            group_by(Index,region, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - Full TS", Difference2="Era 3 - Full TS"))


ratio_upwelling<-rbind(Violin_upwelling%>%filter(period==2)%>%group_by(Index,region, lag, Season)%>%
                       add_column((Violin_upwelling%>%filter(period==3)%>%
                                     group_by(Index,region, lag, Season))$beta-(Violin_upwelling%>%
                                            filter(period==2)%>%
                                            group_by(Index,region, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - 1989:2012", Difference2="Era 3 - Era 2")%>%
  bind_rows(rbind(Violin_upwelling%>%filter(period==3)%>%group_by(Index,region, lag, Season)%>%
                       add_column((Violin_upwelling%>%filter(period==3)%>%
                                     group_by(Index,region, lag, Season))$beta-(Violin_upwelling%>%
                                            filter(period==4)%>%
                                            group_by(Index,region, lag, Season))$beta))%>%
                  rename('beta_diff'=`... - ...`)%>%mutate(Difference="2013:2023 - Full TS", Difference2="Era 3 - Full TS"))


##### Upwelling Violin ####
# plotting functions
q.50 <- function(x) { return(quantile(x, probs=c(0.25,0.75))) }
q.90 <- function(x) { return(quantile(x, probs=c(0.05,0.95))) }

col4 <-pnw_palette(name="Starfish",n=4,type="discrete")
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ratio.up<-ratio_upwelling%>%filter(region!="GoA", Season=="Spring", lag==0)%>%
  mutate(region=ifelse(region=="Northern CC","NCC",
                      ifelse(region=="Southern CC", "SCC", "CCC")))
ratio.up$region <- factor(ratio.up$region,
                levels = c("SCC","CCC","NCC"))
violin.up <-ggplot(ratio.up, aes(x=region, y=beta_diff, fill=region)) +
  theme_bw() +
  scale_fill_manual(values=col4[3:1], name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("Upwelling Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="bottom")
violin.up



unique(ratio_bio$region)
unique(ratio_biological$region)
ratio.bio<-ratio_biological%>%filter(Season=="Spring", lag==0)%>%
         mutate(survey=ifelse(survey=="CALCOFI (SCC)","CALCOFI",
                      ifelse(survey=="RREAS (CCC)", "RREAS",
                        ifelse(survey=="N. Copepod (NCC)", "N. Copepod", 
                               ifelse(survey=="S. Copepod (NCC)","S. Copepod",
                                      ifelse(survey=="Southern Copepod (NCC)","S. Copepod",
                                             ifelse(survey=="Northern Copepod (NCC)","N. Copepod",
                                                    ifelse(survey=="RREAS (SCC)","RREAS",
                                             survey))))))))
ratio.bio$region <- factor(ratio.bio$region,
                levels = c("SCC","CCC","NCC"))
ratio.bio$survey <- factor(ratio.bio$survey,
                levels = c("CALCOFI","RREAS","N. Copepod","S. Copepod"))
violin.bio <-ggplot(ratio.bio%>%filter(Difference=="2013:2023 - 1989:2012"), aes(x=survey, y=beta_diff, fill=region)) +
  theme_bw() +
  scale_fill_manual(values=c(col4[3:1],col4[1]), name="Region")+
  geom_violin(alpha = 0.75, lwd=0.1, scale='width',trim=TRUE) +
  # stat_summary(fun="q.95", colour="black", geom="line", lwd=0.75) +
  stat_summary(fun="q.90", colour="black", geom="line", lwd=0.3) +
  stat_summary(fun="q.50", colour="black", geom="line", lwd=0.6)+
  coord_flip() +
  stat_summary(fun="median", colour="black", size=1, geom="point", pch=21) +
  ggh4x::facet_grid2(Index~Difference2) +
  ylab("Biological Posterior Difference") +
  xlab("") +
  geom_hline(aes(yintercept=0), size=0.3) +
  theme(legend.position="none")



pdf("Output/ViolinV1.pdf", 7,5) 
ggarrange(violin.bio,violin.up,labels = c("A", "B"),  
          font.label = list(size = 12, face="plain"), ncol=2, heights=c(3,0.5))
dev.off()


violin.bio
pdf(file = "Output/Figures/violin.pdf",   # The directory you want to save the file in
    width = 6, # The width of the plot in inches
    height = 6)
violin.up
dev.off()
