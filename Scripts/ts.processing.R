# Timeseries processing and visualization

# PURPOSE: To load and process GOA and EBS stock and climate time series

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(zoo)
library(stats)
library(mgcv)


### LOAD DATA -----------------------------------------------------------------------------------------------------------
  # BSAI and GOA community timeseries
    bsai.ts <- read.csv("./Data/BSAItimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
    #na.omit() # omitting earlier in the ts when all stocks aren't available
    
    goa.ts <- read.csv("./Data/GOAtimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
    #dplyr::select(grep("X", colnames(.), invert = TRUE)) 
    #na.omit() # omitting earlier in the ts when all stocks aren't available
                  

### PROCESS DATA ----------------------------------------------------------------------------------------------------------
  # SST --------------------------
    # BSAI and GOA ERA5 SST
    goa.sst <- read.csv("./Data/sst_GOA.csv") %>%
      rename(Year = year) %>%
      group_by(Year) %>%
      reframe(mean.sst = mean(mean.sst))
    
    ebs.sst <- read.csv("./Data/sst_EBS.csv") %>%
      rename(Year = year) %>%
      group_by(Year) %>%
      reframe(mean.sst = mean(mean.sst))
    
    # bind both datasets
    rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"), 
          goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst
    
  # ADD GF COHORTS AND SALMON OCEAN ENTRY YEAR --------------
    cohorts <- data.frame(Stock = names(c(bsai.ts %>%
                                            dplyr::select(grep(".r0", colnames(.))),
                                          goa.ts %>%
                                            dplyr::select(grep(".r0", colnames(.))))),
                          Lag.Value = c(1, 1, 0, 0, 1, 2, 1, 1, 3, 3, 3, 3, 0, 2, 1, 0, 0, 0, 0, #BSAI
                                        1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 4, 4, 3)) #GOA
    
    ocean.entry <- data.frame(Stock = names(c(bsai.ts %>%
                                                dplyr::select(grep("catch", colnames(.))),
                                              goa.ts %>%
                                                dplyr::select(grep("catch", colnames(.)))))) %>%
                              mutate(Lag.Value = case_when((grepl("sck", Stock) == TRUE) ~ 2,
                                                           (grepl("pnk", Stock) == TRUE) ~ 1,
                                                           (grepl("chm", Stock) == TRUE) ~ 3))   
    
  # BSAI GROUNDFISH SSB ---------------
    # Select  BSAI ssb and log transform
      bsai.ssb <- bsai.ts %>%
        dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
        mutate(bsai.ebs.pol.ssb = as.numeric(bsai.ebs.pol.ssb)) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="SSB") %>%
        mutate(SSB = log(SSB+10), Region = "BSAI") %>%
        rename(Year = YEAR) 
    
  # GOA GROUNDFISH SSB ---------------
    # Select GOA ssb and log transform
      goa.ssb <- goa.ts %>%
        dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="SSB") %>%
        mutate(SSB = log(SSB +10), Region = "GOA") %>%
        rename(Year = YEAR)
    
  # SALMON CATCH ---------------
    # Select catch and log transform
      goa.salmon <- goa.ts %>%
        dplyr::select(YEAR, grep("catch", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="Catch") %>%
        mutate(Catch = log(Catch+10), Region = "GOA") %>%
        rename(Year = YEAR) %>%
        right_join(ocean.entry, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
      bsai.salmon <- bsai.ts %>%
        dplyr::select(YEAR, grep("catch", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="Catch") %>%
        mutate(Catch = log(Catch+10), Region = "BSAI") %>%
        rename(Year = YEAR) %>%
        right_join(ocean.entry, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
      salmon.catch <- rbind(goa.salmon, bsai.salmon)
    
  # CRAB MB ---------------
    # Select BSAI ssb and log transform
      crab.mb <- rbind(bsai.ts %>%
                         dplyr::select(YEAR, grep("mmb", colnames(.))) %>%
                         pivot_longer(., !YEAR, names_to = "Stock", values_to="Value") %>%
                         mutate(Value = log(Value+10), Type = "mmb") %>%
                         rename(Year = YEAR),
                       bsai.ts %>%
                         dplyr::select(YEAR, grep("fmb", colnames(.))) %>%
                         pivot_longer(., !YEAR, names_to = "Stock", values_to="Value") %>%
                         mutate(Value = log(Value+10), Type = "fmb") %>%
                         rename(Year = YEAR)) %>%
        mutate(Stock = case_when((Stock == "bsai.bbrkc.mmb")~ "Bristol Bay red king crab",
                                 (Stock %in% c("bsai.opi.fmb", "bsai.opi.mmb")) ~ "Snow crab",
                                 (Stock %in% c("bsai.tanner.fmb", "bsai.tanner.mmb")) ~ "Tanner crab"))
    
  # GROUNDFISH AND CRAB RECRUITMENT ---------------
    # Select  BSAI r0 and log transform
      bsai.r0 <- bsai.ts %>%
        dplyr::select(YEAR, grep("r0", colnames(.))) %>%
        mutate(bsai.ebs.pol.r0 = as.numeric(bsai.ebs.pol.r0)) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="Recruitment") %>%
        mutate(Recruitment = log(Recruitment +10), Region = "BSAI") %>%
        rename(Year = YEAR) %>%
        right_join(cohorts, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
    # Select  GOA r0 and log transform
      goa.r0 <- goa.ts %>%
        dplyr::select(YEAR, grep("r0", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "Stock", values_to="Recruitment") %>%
        mutate(Recruitment = log(Recruitment +10), Region = "GOA") %>%
        rename(Year = YEAR) %>%
        right_join(cohorts, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
    
  # BSAI GROUNDFISH SSB and SST ---------------
      bsai.ssb %>%
        right_join(., expand.grid(Year = min(bsai.ssb$Year):max(bsai.ssb$Year))) %>% # matching sst and community timeseries by year availability
        right_join(., ebs.sst) %>%
        na.omit() -> bsai.ssb.sst
    
  # GOA GROUNDFISH SSB and SST ---------------
      goa.ssb %>%
        right_join(., expand.grid(Year = min(goa.ssb$Year):max(goa.ssb$Year))) %>%
        right_join(., goa.sst) %>%
        na.omit() -> goa.ssb.sst
      
  # SALMON CATCH AND SST ---------------
      goa.salmon %>%
        right_join(., expand.grid(Year = min(goa.salmon$Year):max(goa.salmon$Year))) %>%
        right_join(., goa.sst %>% rename(Lagged.Year = Year, Lagged.sst = mean.sst), by = c("Lagged.Year")) %>%
        na.omit() -> goa.salmon.sst
      
      bsai.salmon %>%
        right_join(., expand.grid(Year = min(bsai.salmon$Year):max(bsai.salmon$Year))) %>%
        right_join(., ebs.sst %>% rename(Lagged.Year = Year, Lagged.sst = mean.sst), by = c("Lagged.Year")) %>%
        na.omit() -> bsai.salmon.sst
      
      salmon.catch.sst <- rbind(goa.salmon.sst, bsai.salmon.sst)
      
  # CRAB MB AND SST ---------------
      crab.mb %>%
        right_join(., expand.grid(Year = min(crab.mb$Year):max(crab.mb$Year))) %>%
        right_join(., ebs.sst) %>%
        na.omit() -> crab.mb.sst
      
  # RECRUITMENT AND SST ---------------
    # BSAI 
      bsai.r0 %>%
        right_join(., expand.grid(Year = min(bsai.r0$Year):max(bsai.r0$Year))) %>%
        right_join(., ebs.sst %>% rename(Lagged.Year = Year, Lagged.sst = mean.sst), by = c("Lagged.Year")) %>%
        na.omit() -> bsai.r0.sst
      
    # GOA
      goa.r0 %>%
        right_join(., expand.grid(Year = min(goa.r0$Year):max(goa.r0$Year))) %>%
        right_join(., goa.sst %>% rename(Lagged.Year = Year, Lagged.sst = mean.sst), by = c("Lagged.Year")) %>%
        na.omit() -> goa.r0.sst
      
    
    
      
### FIND YEAR RANGES WHERE ALL TS ARE PRESENT ---------------------------------------------------
  bsai.ssb %>%
      na.omit() -> tt
    
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1)) -> hh #1992:2021
  
  goa.ssb %>%
      na.omit() -> tt
  
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1)) -> hh #1982:2021

  bsai.salmon %>%
    na.omit() -> tt
  
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1))-> hh #1893:2022
  
  
  goa.salmon %>%
      na.omit() %>%
      filter(Year > 1980)-> tt
      
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1))-> hh #1995:2021
      
  crab.mb %>%
     filter(Type == "mmb") %>%
     na.omit() -> tt
  
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1))-> hh #1985:2021
      
  crab.mb %>%
      filter(Type == "fmb") %>%
      na.omit() -> tt
      
      hist(tt$Year, breaks = seq(min(tt$Year), max(tt$Year), by = 1))-> hh #1982:2021
      
  years <- 1995:2021
  
    