# Timeseries processing and visualization

# PURPOSE: To load and process GOA and EBS TS and climate time series

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(zoo)
library(stats)
library(mgcv)
library(MuMIn)
library(roll)
library(TTR)

### LOAD DATA -----------------------------------------------------------------------------------------------------------
  # BSAI and GOA community timeseries
    bsai.ts <- read.csv("./Data/BSAItimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
    #na.omit() # omitting earlier in the ts when all TS aren't available
    
    goa.ts <- read.csv("./Data/GOAtimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
    #dplyr::select(grep("X", colnames(.), invert = TRUE)) 
    #na.omit() # omitting earlier in the ts when all TSs aren't available
                  

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
    cohorts <- data.frame(TS = names(c(bsai.ts %>%
                                            dplyr::select(grep(".r0", colnames(.))),
                                          goa.ts %>%
                                            dplyr::select(grep(".r0", colnames(.))))),
                          Lag.Value = c(1, 1, 0, 0, 1, 2, 1, 1, 3, 3, 3, 3, 0, 2, 1, 0, 0, 0, 0, #BSAI
                                        1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 4, 4, 3)) #GOA
    
    ocean.entry <- data.frame(TS = names(c(bsai.ts %>%
                                                dplyr::select(grep("catch", colnames(.))),
                                              goa.ts %>%
                                                dplyr::select(grep("catch", colnames(.)))))) %>%
                              mutate(Lag.Value = case_when((grepl("sck", TS) == TRUE) ~ 2,
                                                           (grepl("pnk", TS) == TRUE) ~ 1,
                                                           (grepl("chm", TS) == TRUE) ~ 3))   
    
  # BSAI GROUNDFISH SSB ---------------
    # Select  BSAI ssb and log transform
      bsai.ssb <- bsai.ts %>%
        dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
        mutate(bsai.ebs.pol.ssb = as.numeric(bsai.ebs.pol.ssb)) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="SSB") %>%
        mutate(log.SSB = log(SSB+10), SSB = SSB, Region = "BSAI") %>%
        rename(Year = YEAR) 
    
  # GOA GROUNDFISH SSB ---------------
    # Select GOA ssb and log transform
      goa.ssb <- goa.ts %>%
        dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="SSB") %>%
        mutate(log.SSB = log(SSB+10), SSB = SSB, Region = "GOA") %>%
        rename(Year = YEAR)
    
  # SALMON CATCH ---------------
    # Select catch and log transform
      goa.salmon <- goa.ts %>%
        dplyr::select(YEAR, grep("catch", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="Catch") %>%
        mutate(log.catch = log(Catch+10), Catch = Catch, Region = "GOA") %>%
        rename(Year = YEAR) %>%
        right_join(ocean.entry, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
      bsai.salmon <- bsai.ts %>%
        dplyr::select(YEAR, grep("catch", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="Catch") %>%
        mutate(log.catch = log(Catch+10), Catch = Catch, Region = "BSAI") %>%
        rename(Year = YEAR) %>%
        right_join(ocean.entry, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
      salmon.catch <- rbind(goa.salmon, bsai.salmon)
    
  # CRAB MB ---------------
    # Select BSAI ssb and log transform
      crab.mb <- rbind(bsai.ts %>%
                         dplyr::select(YEAR, grep("mmb", colnames(.))) %>%
                         pivot_longer(., !YEAR, names_to = "TS", values_to="Value") %>%
                         mutate(log.value = log(Value+10), Value = Value, Type = "mmb") %>%
                         rename(Year = YEAR),
                       bsai.ts %>%
                         dplyr::select(YEAR, grep("fmb", colnames(.))) %>%
                         pivot_longer(., !YEAR, names_to = "TS", values_to="Value") %>%
                         mutate(log.value = log(Value+10), Value = Value, Type = "fmb") %>%
                         rename(Year = YEAR)) %>%
        mutate(TS = case_when((TS == "bsai.bbrkc.mmb")~ "Bristol Bay red king crab",
                                 (TS %in% c("bsai.opi.fmb", "bsai.opi.mmb")) ~ "Snow crab",
                                 (TS %in% c("bsai.tanner.fmb", "bsai.tanner.mmb")) ~ "Tanner crab"))
    
  # GROUNDFISH AND CRAB RECRUITMENT ---------------
    # Select  BSAI r0 and log transform
      bsai.r0 <- bsai.ts %>%
        dplyr::select(YEAR, grep("r0", colnames(.))) %>%
        mutate(bsai.ebs.pol.r0 = as.numeric(bsai.ebs.pol.r0)) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="Recruitment") %>%
        mutate(log.recruitment = log(Recruitment +10), Recruitment = Recruitment, Region = "BSAI") %>%
        rename(Year = YEAR) %>%
        right_join(cohorts, .) %>%
        mutate(Lagged.Year = Year - Lag.Value)
      
    # Select  GOA r0 and log transform
      goa.r0 <- goa.ts %>%
        dplyr::select(YEAR, grep("r0", colnames(.))) %>%
        pivot_longer(., !YEAR, names_to = "TS", values_to="Recruitment") %>%
        mutate(log.recruitment = log(Recruitment +10), Recruitment = Recruitment, Region = "GOA") %>%
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
      
    
    
      
### MAKE LABELS -------------------------------------------------------------------
      # Specify labels
      ssb.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                         "Atka mackerel", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                         "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                         "Pacific ocean perch", "Sablefish", "Skate", "Greenland turbot", "Yellowfin sole")
      
      names(ssb.labs.bsai) <- c("bsai.ai.cod.ssb", "bsai.ai.pol.ssb", "bsai.apl.ssb", "bsai.atf.ssb", "bsai.atk.ssb",
                                "bsai.brr.ssb", "bsai.ebs.cod.ssb", "bsai.ebs.pol.ssb", "bsai.kam.ssb", "bsai.nrf.ssb",
                                "bsai.nrs.ssb", "bsai.pop.ssb", "bsai.sab.ssb", "bsai.ska.ssb", "bsai.turb.ssb", 
                                "bsai.yfs.ssb")
      
      # Specify labels
      ssb.labs.goa <- c("Arrowtooth flounder", "Cod", "Rougheye/blackspotted rockfish", "Dover sole", "Dusky rockfish",
                        "Flathead sole", "Northern rockfish", "Northern rock sole", "Pollock", "Pacific ocean perch",
                        "Rex sole", "Sablefish", "Southern rock sole")
      
      names(ssb.labs.goa) <- c("goa.atf.ssb", "goa.cod.ssb", "goa.dbr.ssb", "goa.dov.ssb", "goa.drf.ssb", "goa.fhs.ssb",
                               "goa.nrf.ssb", "goa.nrs.ssb", "goa.pol.ssb", "goa.pop.ssb", "goa.rex.ssb", "goa.sab.ssb",
                               "goa.srs.ssb")
      
      # Specify labels
      salm.labs <- c("Bristol Bay sockeye", "Chignik chum", "Chignik pink", "Chignik sockeye", "Cook Inlet chum",
                     "Cook Inlet pink", "Cook Inlet sockeye", "Kodiak chum", "Kodiak pink", "Kodiak sockeye",
                     "Peninsula chum", "Peninsula pink", "Peninsula sockeye", "PWS chum",
                     "PWS pink", "PWS sockeye", "Southeast chum", "Southeast pink",
                     "Southeast sockeye")
      
      names(salm.labs) <- c("bsai.sck.catch", "chig.chm.catch", "chig.pnk.catch", "chig.sck.catch",
                            "cook.chm.catch", "cook.pnk.catch", "cook.sck.catch", "kod.chm.catch",
                            "kod.pnk.catch", "kod.sck.catch", "pen.chm.catch", "pen.pnk.catch",
                            "pen.sck.catch", "pws.chm.catch", "pws.pnk.catch", "pws.sck.catch", "se.chm.catch",
                            "se.pnk.catch", "se.sck.catch")
      
      # Specify BSAI labels
      r0.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                        "Atka mackerel", "Bristol Bay red king crab", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                        "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                        "Snow crab", "Pacific ocean perch", "Sablefish", "Skate", 
                        "Tanner crab", "Greenland turbot", "Yellowfin sole")
      
      names(r0.labs.bsai) <- c("bsai.ai.cod.r0", "bsai.ai.pol.r0", "bsai.apl.r0", "bsai.atf.r0", "bsai.atk.r0",
                               "bsai.bbrkc.r0", "bsai.brr.r0", "bsai.ebs.cod.r0", "bsai.ebs.pol.r0", "bsai.kam.r0", "bsai.nrf.r0",
                               "bsai.nrs.r0", "bsai.opi.r0", "bsai.pop.r0", "bsai.sab.r0", "bsai.ska.r0", 
                               "bsai.tanner.r0", "bsai.turb.r0",  "bsai.yfs.r0")
      
      # Specify GOA labels
      r0.labs.goa <- c("Arrowtooth flounder", "Cod", "Rougheye/blackspotted rockfish", "Dover sole", "Dusky rockfish",
                       "Flathead sole", "Northern rockfish", "Northern rock sole", "Pollock", "Pacific ocean perch",
                       "Rex sole", "Sablefish", "Southern rock sole")
      
      names(r0.labs.goa) <- c("goa.atf.r0", "goa.cod.r0", "goa.dbr.r0", "goa.dov.r0", "goa.drf.r0", "goa.fhs.r0",
                              "goa.nrf.r0", "goa.nrs.r0", "goa.pol.r0", "goa.pop.r0", "goa.rex.r0", "goa.sab.r0",
                              "goa.srs.r0")
      
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
      
  years <- 1995:2021 # years where all timeseries are present
  
    