# Data

## Project contacts
Megan Feddern, mfeddern@alaska.edu
Eric Ward, eric.ward@noaa.gov
Mary Hunsicker, mary.hunsicker@noaa.gov
Will Satterthwaite, will.satterthwaite@noaa.gov

### DATA-SPECIFIC INFORMATION FOR: all_juvenile_indices.rds

Summary: Juvenile indices of groundfish abundance from 2003 - 2022 (missing 2020 due to the pandemic). Indices are derived from sdmTMB which was run by Nick Bond data is for CCE. 

1. Number of variables: 7
2. Number of cases/rows: 211
3. Variable List: 
    year: year abundance was measured
    est: estimate of abundance
    lwr: lower bound of uncertainty for est
    upr: upper bound of uncertainty for est
    log_est: log transformed estimate of abundance
    se: standard error of log_est
    common_name: common name of species with estimated abundance
4. Missing data codes: 
    NA
5. Abbreviations used: 
    SSB; Spawning Stock Biomass
6. Other relevant information:
    Data is from: 
    Species list: arrowtooth flounder, darkblotched rockfish, dover sole, lingcod, longspine thornyhead, Pacific grenadier, Pacific hake, Pacific sanddab, sablefish, shortspine thornyhead, splitnose rockfish
    
### DATA-SPECIFIC INFORMATION FOR: biologydata_south.central_2023-All

Summary: standardized abundance indices for central (RREAS) and southern (CALCOFI core area) updated to 2022 and standardized according to the methods described in Hunsicker et al. 2022 "Tracking and forecasting community responses to climate perturbations in the California Current Ecosystem". All other columns list a species name, each cell contains an index of abundance for species listed in column heading for a given year
  
  1. Number of variables: 45
  2. Number of cases/rows: 73
  3. Variable List: 
      year: year abundance was measured
      ZALOPHUS.PUPCT: Average number of California sea lion pups at San Miguel Island for cohorts from 1997-present. (Note: In 2020, data collection was limited to estimates from aerial surveys due to COVID-19 restrictions to field operations. M. Ball [Wildlands Conservation Science] conducted the aerial surveys and E. Jaime [AFSC] interpreted images to derive counts.) Source: AFSC/NMML (http://www.afsc.noaa.gov/nmml/species/species_cal.php) and summarized in CCIEA dashboard
      ZALOPHUS.PUPWT: Predicted average pup weights for female California sea lion pups born at San Miguel Island, California. Pups are weighed in September or October each year and weights are adjusted using a mixed effects model to a 1 October weighing date.  Source: AFSC/NMML (http://www.afsc.noaa.gov/nmml/species/species_cal.php) and summarized in CCIEA dashboard
      ZALOPHUS.PUPGROWTH: predicted average daily growth rate of female California sea lion pups at San Miguel Island between 4-7 months of age for cohorts 1997 - present. Source: AFSC/NMML (http://www.afsc.noaa.gov/nmml/species/species_cal.php) and summarized in CCIEA dashboard
      SBRD.ASSP.PROD: Central CC - SEFI Ashy Storm Petrel productivity anomaly. 
      SBRD.BRCO.PROD:  Central CC - SEFI Brandt's cormorant productivity anomaly. 
      SBRD.CAAU.PROD: Central CC  - SEFI Cassin's auklet productivity anomaly
      SBRD.COMU.PROD: Ce CC - SEFI Common murre productivity anomaly
      SBRD.PECO.PROD: Central CC - SEFI Pelagic Cormorant productivity anomaly. 
      SBRD.PIGO.PROD: Ce CC - SEFI Pigeon guillemot productivity anomaly
      SBRD.RHAU.PROD: Ce CC - SEFI Rhinoceros auklet productivity anomaly
      SBRD.WEGO.PROD: Central CC - SEFI Western Gull productivity anomaly. 
      rreas.adu.anchovy: 1990 - 2022 adult Northern anchovy from RREAS survey
      rreas.adu.sardine: 1990 - 2022 adult Pacific sarding from RREAS survey
      rreas.all.yoy.rockfish: 1983 - 2022 Rockfish young of year for all species RREAS survey
      rreas.krill: 1990 - 2022 krill from  RREAS survey
      rreas.market.squid: 1990 - 2022 Market squid from RREAS survey
      rreas.myctophids: 1990 - 2022 myctophids (Lanternfishes) from RREAS survey
      rreas.octopus: 1990 - 2022 octupus from RREAS survey
      rreas.yoy.anchovy: 1990 - 2022 young of year Northern anchovy from RREAS survey
      rreas.yoy.bocaccio: 1984 - 2022 young of year Bocaccio rockfish from RREAS survey
      rreas.yoy.chili: 1984 - 2022 young of year chilipepper rockfish from RREAS survey
      rreas.yoy.hake: 1983 - 2022 young of year Pacific hake from RREAS survey
      rreas.yoy.halfbanded: 1984 - 2022 young of year halfbanded rockfish from RREAS survey
      rreas.yoy.lingcod: 1984 - 2022 young of year lingcod from RREAS survey
      rreas.yoy.pacdabs: 1987 - 2022 young of year Pacific sandabs from RREAS survey
      rreas.yoy.sardine: 1995 - 2022 young of year Pacific sardine from RREAS survey Note: only 9 years of data - multiple missing years
      rreas.yoy.shortbelly: 1983 - 2022 young of year shortbelly rockfish from RREAS survey
      rreas.yoy.speckdabs: 1987 - 2022 young of year speckled sanddabs from RREAS survey
      rreas.yoy.widow: 1983 - 2022 young of year widow rockfish from RREAS survey
      rreas.yoy.ytail: 1983 - 2022 young of year yellowtail rockkfish from RREAS survey
      calcofi.bathylagoides.wesethi: 1951 - 2022 snubnose smelt from CalCOFI survey
      calcofi.bathylagus.pacificus: 1951 - 2022 slender blacksmelt from CalCOFI survey
      calcofi.ceratoscopelus.townsendi: 1951 - 2022 dogtooth lampfish from CalCOFI survey
      calcofi.engraulis.mordax: 1951 - 2022 northern anchovy from CalCOFI survey
      calcofi.leuroglossus.stilbius: 1951 - 2022 California smoothtongue from CalCOFI survey
      calcofi.lipolagus.ochotensis: 1951 - 2022 eared blacksmelt from CalCOFI survey
      calcofi.merluccius.productus: 1951 - 2022 north Pacific hake from CalCOFI survey
      calcofi.protomyctophum.crockeri: 1951 - 2022 California flashlightfish from CalCOFI survey
      calcofi.sardinops.sagax: 1951 - 2022 sardinops from CalCOFI survey
      calcofi.stenobrachius.leucopsarus: 1951 - 2022 northern lampfish from CalCOFI survey
      calcofi.symbolophorus.californiensis: 1951 - 2022 Bigfin lanternfish from CalCOFI survey
      calcofi.tarletonbeania.crenularis: 1951 - 2022  Blue lanternfish from CalCOFI survey
      calcofi.triphoturus.mexicanus: 1951 - 2022 Mexican lampfish from CalCOFI survey
      calcofi.vinciguerria1: 1951 - 2022 lightfishes from CalCOFI survey
  4. Missing data codes: 
      NA
  5. Abbreviations used: 
  6. Other relevant information:
      Data is from: CalCOFI and RREAs; summarized and standardized by Mary Hunsicker, can be accessed publically from the ERDAPP server (https://coastwatch.pfeg.noaa.gov/erddap/index.html)
      Sea bird data: Data from Point Blue Conservation Science collected on Southeast Farallon Island in collaboration with the Farallon Islands National Wildlife Refuge (USFWS); contact Dr. Jaime Jahncke (jjahncke@pointblue.org) before citing or distributing these data. Productivity anomaly is the annual mean number of chicks fledged per breeding pair per species minus the long term mean, which is calculated by averaging all of the annual means prior to the most recent year (for data from 1986 to 2018, the long term mean is calculated including data from 1986-2017. (https://www.integratedecosystemassessment.noaa.gov/regions/california-current)

