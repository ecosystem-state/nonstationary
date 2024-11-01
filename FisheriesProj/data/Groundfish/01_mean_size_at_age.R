
library(tidyverse)
library(ggplot2)
library(sp)
library(sf)
# library(rgdal)
# remotes::install_github("pfmc-assessments/nwfscSurvey", 'test-sample-type', force=TRUE)
library(nwfscSurvey)

# pull all biological samples
fish <- PullBio.fn(SurveyName = "NWFSC.Combo")
names(fish) <- tolower(names(fish))

sub <- dplyr::filter(fish, age <=2, sex=="F")

# aggregate by species and year
min_sample <- 5
min_years <- 10
summaries <- dplyr::group_by(sub, common_name, year, age) %>%
  dplyr::summarise(n = n(),
                   scientific_name = scientific_name[1],
                   mean_length = mean(length_cm,na.rm=T),
                   mean_weight = mean(weight,na.rm=T)) %>%
  dplyr::filter(n >= min_sample) %>%
  dplyr::group_by(common_name, age) %>%
  dplyr::mutate(nyr = length(unique(year))) %>%
  dplyr::filter(nyr > min_years) %>%
  dplyr::select(-nyr)

saveRDS(summaries, "wcbts_mean_size_at_age.rds")

