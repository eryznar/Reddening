#Timeseries processing and visualization

#PURPOSE: To load, process, and visualize GOA and EBS stock and climate timeseries

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)

### LOAD DATA -----------------------------------------------------------------------------------------------------------
bsai.ts <- read.csv("./Data/BSAItimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
              #na.omit() # omitting earlier in the ts when all stocks aren't available
goa.ts <- read.csv("./Data/GOAtimeseries.csv", stringsAsFactors=FALSE, fileEncoding="latin1") 
  dplyr::select(grep("X", colnames(.), invert = TRUE)) 
  #na.omit() # omitting earlier in the ts when all stocks aren't available
    
goa.sst <- read.csv("./Data/sst_GOA.csv") %>%
  rename(Year = year) %>%
  group_by(Year) %>%
  reframe(mean.sst = mean(mean.sst))

ebs.sst <- read.csv("./Data/sst_EBS.csv") %>%
  rename(Year = year) %>%
  group_by(Year) %>%
  reframe(mean.sst = mean(mean.sst))


### PROCESS/PLOT RAW TS DATA ----------------------------------------------------------------------------------------------
  # BSAI GROUNDFISH --------
    # Select  BSAI ssb and log transform
    bsai.ssb <- bsai.ts %>%
      dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
      mutate(bsai.ebs.pol.ssb = as.numeric(bsai.ebs.pol.ssb)) %>%
      pivot_longer(., !YEAR, names_to = "Stock", values_to="SSB") %>%
      mutate(SSB = log(SSB), Region = "BSAI") %>%
      rename(Year = YEAR)
    
    # Specify labels
    ssb.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                  "Atka mackerel", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                  "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                  "Pacific ocean perch", "Sablefish", "Skate", "Greenland turbot", "Yellowfin sole")
    
    names(ssb.labs.bsai) <- c("bsai.ai.cod.ssb", "bsai.ai.pol.ssb", "bsai.apl.ssb", "bsai.atf.ssb", "bsai.atk.ssb",
                         "bsai.brr.ssb", "bsai.ebs.cod.ssb", "bsai.ebs.pol.ssb", "bsai.kam.ssb", "bsai.nrf.ssb",
                         "bsai.nrs.ssb", "bsai.pop.ssb", "bsai.sab.ssb", "bsai.ska.ssb", "bsai.turb.ssb", 
                         "bsai.yfs.ssb")
    # Plot
    ggplot(bsai.ssb %>% filter(Year > 1952), aes(x = Year, y = SSB))+
      geom_line(linewidth = 1)+
      facet_wrap(~Stock, scales = "free_y", labeller = labeller(Stock = ssb.labs.bsai))+
      theme_bw()+
      ggtitle("BSAI groundfish spawning stock biomass")+
      ylab("log(kilotons)")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10))-> bsai.ssb.plot
    
    ggsave(plot = bsai.ssb.plot, "./Figures/bsai.ssb.plot.png", width = 12, height = 8, units = "in")
    
    # Select  BSAI r0 and log transform
    bsai.r0 <- bsai.ts %>%
      dplyr::select(YEAR, grep("r0", colnames(.))) %>%
      mutate(bsai.ebs.pol.r0 = as.numeric(bsai.ebs.pol.r0)) %>%
      pivot_longer(., !YEAR, names_to = "Stock", values_to="Recruitment") %>%
      mutate(Recruitment = log(Recruitment), Region = "BSAI") %>%
      rename(Year = YEAR)
    
    # Specify labels
    r0.labs.bsai <- c("Aleutian Islands cod", "Aleutian Islands pollock", "Alaska plaice", "Arrowtooth flounder",
                       "Atka mackerel", "Bristol Bay red king crab", "Blackspotted/rougheye rockfish", "Eastern Bering Sea cod",
                       "Eastern Bering Sea pollock", "Kamchatka flounder", "Northern rockfish", "Northern rock sole",
                       "Snow crab", "Pacific ocean perch", "Sablefish", "Skate", 
                      "Tanner crab", "Greenland turbot", "Yellowfin sole")
    
    names(r0.labs.bsai) <- c("bsai.ai.cod.r0", "bsai.ai.pol.r0", "bsai.apl.r0", "bsai.atf.r0", "bsai.atk.r0",
                              "bsai.bbrkc.r0", "bsai.brr.r0", "bsai.ebs.cod.r0", "bsai.ebs.pol.r0", "bsai.kam.r0", "bsai.nrf.r0",
                              "bsai.nrs.r0", "bsai.opi.r0", "bsai.pop.r0", "bsai.sab.r0", "bsai.ska.r0", 
                             "bsai.tanner.r0", "bsai.turb.r0",  "bsai.yfs.r0")
    # Plot
    ggplot(bsai.r0 %>% filter(Year > 1950), aes(x = Year, y = Recruitment))+
      geom_line(linewidth = 0.75)+
      facet_wrap(~Stock, scales = "free_y", labeller = labeller(Stock = r0.labs.bsai))+
      theme_bw()+
      ggtitle("BSAI groundfish and crab recruitment")+
      ylab("log(millions of recruits)") + 
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10))-> bsai.r0.plot
    
    ggsave(plot = bsai.r0.plot, "./Figures/bsai.r0.plot.png", width = 12, height = 8, units = "in")
  
  # GOA GROUNDFISH -----------
    # Select GOA ssb and log transform
    goa.ssb <- goa.ts %>%
      dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
      pivot_longer(., !YEAR, names_to = "Stock", values_to="SSB") %>%
      mutate(SSB = log(SSB), Region = "GOA") %>%
      rename(Year = YEAR)
    
    # Specify labels
    ssb.labs.goa <- c("Arrowtooth flounder", "Cod", "Dusky/blackspotted rockfish", "Dover sole", "Dusky rockfish",
                  "Flathead sole", "Northern rockfish", "Northern rock sole", "Pollock", "Pacific ocean perch",
                  "Rex sole", "Sablefish", "Southern rock sole")
    
    names(ssb.labs.goa) <- c("goa.atf.ssb", "goa.cod.ssb", "goa.dbr.ssb", "goa.dov.ssb", "goa.drf.ssb", "goa.fhs.ssb",
                         "goa.nrf.ssb", "goa.nrs.ssb", "goa.pol.ssb", "goa.pop.ssb", "goa.rex.ssb", "goa.sab.ssb",
                         "goa.srs.ssb")
    
    # Plot
    ggplot(goa.ssb %>% filter(Year > 1959), aes(x = Year, y = SSB))+
      geom_line(linewidth = 1)+
      facet_wrap(~Stock, scales = "free_y", labeller = labeller(Stock = ssb.labs.goa))+
      theme_bw()+
      ggtitle("GOA groundfish spawning stock biomass")+
      ylab("log(kilotons)")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10)) -> goa.ssb.plot
    
    ggsave(plot = goa.ssb.plot, "./Figures/goa.ssb.plot.png", width = 12, height = 8, units = "in")
    
    
  # SALMON ----------------
    # Select catch and log transform
    goa.salmon <- goa.ts %>%
      dplyr::select(YEAR, grep("catch", colnames(.))) %>%
      pivot_longer(., !YEAR, names_to = "Stock", values_to="Catch") %>%
      mutate(Catch = log(Catch), Region = "GOA") %>%
      rename(Year = YEAR)
    
    bsai.salmon <- bsai.ts %>%
      dplyr::select(YEAR, grep("catch", colnames(.))) %>%
      pivot_longer(., !YEAR, names_to = "Stock", values_to="Catch") %>%
      mutate(Catch = log(Catch), Region = "BSAI") %>%
      rename(Year = YEAR)
    
    salmon.catch <- rbind(goa.salmon, bsai.salmon)
    
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
    
    # Plot
    ggplot(salmon.catch %>% filter(Year > 1959), aes(x = Year, y = Catch))+
      geom_line(linewidth = 0.75)+
      facet_wrap(~Stock, scales = "free_y", labeller = labeller(Stock = salm.labs))+
      theme_bw()+
      ggtitle("BSAI and GOA salmon catch")+
      ylab("log(number of fish)")+
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12)) -> salmon.catch.plot
    
    ggsave(plot = salmon.catch.plot, "./Figures/salmon.catch.plot.png", width = 12, height = 8, units = "in")
    
    
  # CRAB ----------------
    # Select  BSAI ssb and log transform
    crab.mb <- rbind(bsai.ts %>%
                      dplyr::select(YEAR, grep("mmb", colnames(.))) %>%
                      pivot_longer(., !YEAR, names_to = "Stock", values_to="Value") %>%
                      mutate(Value = log(Value), Type = "mmb") %>%
                      rename(Year = YEAR),
                    bsai.ts %>%
                      dplyr::select(YEAR, grep("fmb", colnames(.))) %>%
                      pivot_longer(., !YEAR, names_to = "Stock", values_to="Value") %>%
                      mutate(Value = log(Value), Type = "fmb") %>%
                      rename(Year = YEAR)) %>%
              mutate(Stock = case_when((Stock == "bsai.bbrkc.mmb")~ "Bristol Bay red king crab",
                                       (Stock %in% c("bsai.opi.fmb", "bsai.opi.mmb")) ~ "Snow crab",
                                       (Stock %in% c("bsai.tanner.fmb", "bsai.tanner.mmb")) ~ "Tanner crab"))
    
    
    # Plot
    ggplot(crab.mb %>% filter(Year > 1953), aes(x = Year, y = Value, linetype = Type))+
      geom_line(linewidth = 0.75)+
      facet_wrap(~Stock, scales = "free_y", nrow = 3)+
      theme_bw()+
      scale_linetype_manual(name = "", values = c("solid", "dashed"), labels = c("Female", "Male"))+
      ggtitle("BSAI crab mature biomass")+
      ylab("log(kilotons)")+
      theme(legend.position = "bottom",
            axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14),
                  strip.text = element_text(size = 10),
            legend.text = element_text(size = 12)) -> crab.mb.plot
    
    ggsave(plot = crab.mb.plot, "./Figures/crab.mb.plot.png", width = 8, height = 10, units = "in")
    
   
 # SST ---------------
    # bind both datasets
    rbind(ebs.sst %>% mutate(region = "Eastern Bering Sea"), 
          goa.sst %>% mutate(region = "Gulf of Alaska")) -> sst
    
    # plot
    ggplot(sst %>% filter(Year < 2024), aes(x = Year, y = mean.sst))+
      facet_wrap(~region, scales = "free_y", nrow = 2)+
      geom_line(linewidth = 0.75)+
      geom_point(size = 1.5)+
      scale_x_continuous(breaks = seq(min(sst$Year), max(sst$Year), by = 10))+
      theme_bw()+
      ggtitle("ERA5 SST")+
      ylab("°C")+
      theme(legend.position = "bottom",
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            strip.text = element_text(size = 10),
            legend.text = element_text(size = 12)) -> sst.plot
    
    ggsave(plot= sst.plot, "./Figures/ERA5sst.png", width = 6, height = 8, units = "in")
    
### PLOT SSB/CATCH DATA WITH SST ------------
  bsai.ssb %>%
    right_join(., expand.grid(Year = min(bsai.ssb$Year):max(bsai.ssb$Year))) %>%
    right_join(., ebs.sst) %>%
    na.omit() -> bsai.ssb.sst
  
    
  
  # Plot
  ggplot(bsai.ssb.sst %>% filter(Year > 1952), aes(x = mean.sst, y = SSB))+
    geom_line(linewidth = 1)+
    facet_wrap(~Stock, scales = "free_y", labeller = labeller(Stock = ssb.labs.bsai))+
    theme_bw()+
    ggtitle("BSAI groundfish spawning stock biomass")+
    ylab("log(kilotons)")+
    xlab("SST")+
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 10))-> bsai.ssb.sst.plot
  
