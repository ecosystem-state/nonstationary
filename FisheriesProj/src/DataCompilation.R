library(dplyr)
library(reshape2)
library(bayesdfa)
library(MCMCvis)
library(ggplot2)
library(stringr)
library(ggpubr)
library(tidyverse)
library(dplyr)

dat<- readRDS("data/all_juvenile_indices.rds")
dat<- read.csv("data/salmon/MatRate_Survival_up to CY 2021_new.csv")

Indicator_ID<-data.frame(stock_code=unique(dat$stock_code))
Indicator_ID%>%mutate(stock_name=ifelse(stock_code=='ATN', "Atnarko", 
                                 ifelse(stock_code=='BQR', "Big Qualicum River Fall",
                                 ifelse(stock_code=='CHI', "Chilliwack River Fall",
                                 ifelse(stock_code=='COW', "Cowichan River Fall",
                                 ifelse(stock_code=='CWF', "Cowlitz Fall Tule",
                                 ifelse(stock_code=='ELK', "Elk River",
                                 ifelse(stock_code=='ELW', "Elwha River",
                                 ifelse(stock_code=='GAD', "George Adams Fall Fingerling",
                                 ifelse(stock_code=='GRN', "Green River Fingerling",
                                 ifelse(stock_code=='HAN', "Hanford Wild Brights",
                                 ifelse(stock_code=='HAR', "Harrison River",
                                 ifelse(stock_code=='HOK', "Hoko Fall Fingerling",
                                 ifelse(stock_code=='KLM', "Kitsumkalum River Summer",
                                 ifelse(stock_code=='KLY', "Kitsumkalum Yearling",
                                 ifelse(stock_code=='LRH', "Columbia Lower River Hatchery",
                                 ifelse(stock_code=='LRW', "Lewis River Wild",
                                 ifelse(stock_code=='LYF', "Lyons Ferry Fingerling",
                                 ifelse(stock_code=='LYY', "Lyons Ferry Yearling", 
                                 ifelse(stock_code=='MSH', "Middle Shuswap",
                                 ifelse(stock_code=='NIC', "Nicola",
                                 ifelse(stock_code=='NIS', "Nisqually Fall Fingerling",
                                 ifelse(stock_code=='NSA', "Northern Southeast Alaska",
                                 ifelse(stock_code=='NSF', "Nooksack Spring Fingerling",
                                 ifelse(stock_code=='PHI', "Phillips River Fall",
                                 ifelse(stock_code=='PPS', "Puntledge River Summer",
                                 ifelse(stock_code=='QUE', "Queets Fall Fingerling",
                                 ifelse(stock_code=='QUI', "Quinsam River Fall ",
                                 ifelse(stock_code=='RBT', "Robertson Creek Fall",
                                 ifelse(stock_code=='SAM', "Samish Fall Fingerling",
                                 ifelse(stock_code=='SHU', "Lower Shuswap River Summer",
                                 ifelse(stock_code=='SKF', "Skagit Spring Fingerling",
                                 ifelse(stock_code=='SKY', "Skykomish Fall Fingerling",
                                 ifelse(stock_code=='SMK', "Similkameen Summer Yearling",
                                 ifelse(stock_code=='SOO', "Tsoo-Yess Fall Fingerling",
                                 ifelse(stock_code=='SPR', "Spring Creek Tule",
                                 ifelse(stock_code=='SPS', "South Puget Sound Fall Fingerling",
                                 ifelse(stock_code=='SRH', "Salmon River",
                                 ifelse(stock_code=='SSA', "Southern SEAK Spring",
                                 ifelse(stock_code=='SSF', "Skagit Summer Fingerling",
                                 ifelse(stock_code=='STL', "Stillaguamish Fall Fingerling² ",
                                 ifelse(stock_code=='SUM', "Columbia River Summers",
                                 ifelse(stock_code=='URB', "Columbia Upriver Bright",
                                 ifelse(stock_code=='WSH', "Willamette Spring",1))))))))))))))))))))))))))))))))))))))))))))


Indicator_ID%>%mutate(area=ifelse(stock_code=='NSA'|stock_code=='SSA'|stock_code=='CHK'|stock_code=='UNU', "Southeast Alaska", 
                                 ifelse(stock_code=='TST', "Transboundary Rivers",
                                 ifelse(stock_code=='ATN'|stock_code=='KLM'|stock_code=='KLY', "North/Central BC",
                                 ifelse(stock_code=='RBT', "WCVI",
                                 ifelse(stock_code=='QUI'|stock_code=='PHI'|stock_code=='PPS'|stock_code=='BQR'|stock_code=='COW', "Strait of Georgia",
                                 ifelse(stock_code=='HAR'|stock_code=='CHI'|stock_code=='CKO'|stock_code=='NIC'|stock_code=='SHU'|stock_code=='MSH', "Fraser River",
                                 ifelse(stock_code=='NSF'|stock_code=='SAM'|stock_code=='SSF'|stock_code=='SKF', "North Puget Sound",
                                 ifelse(stock_code=='STL'|stock_code=='SKY', "Central Puget Sound",
                                 ifelse(stock_code=='NIS'|stock_code=='SPS'|stock_code=='GRN', "South Puget Sound",
                                 ifelse(stock_code=='GAD', "Hood Canal",
                                 ifelse(stock_code=='ELW', "Juan de Fuca",
                                 ifelse(stock_code=='HOK'|stock_code=='QUE'|stock_code=='SOO', "North Washington Coast",
                                 ifelse(stock_code=='LRH'|stock_code=='CWF'|stock_code=='LRW'|stock_code=='WSH'|stock_code=='SPR', "Lower Columbia River",
                                 ifelse(stock_code=='HAN'|stock_code=='SMK'|stock_code=='SUM'|stock_code=='URB', "Upper Columbia River",
                                 ifelse(stock_code=='LYF'|stock_code=='LYY', "Snake River",
                                 ifelse(stock_code=='SRH', "North Oregon Coast",
                                 ifelse(stock_code=='ELK', "Mid Oregon Coast",1))))))))))))))))))
