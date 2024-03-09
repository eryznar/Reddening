#Timeseries processing and visualization

#PURPOSE: To load, process, and visualize GOA and EBS stock and climate timeseries

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)

### LOAD DATA -----------------------------------------------------------------------------------------------------------
bsai.ts <- read.csv("./Data/BSAItimeseries.csv")
goa.ts <- read.csv("./Data/GOAtimeseries.csv") %>%
  dplyr::select(grep("X", colnames(.), invert = TRUE))

### PROCESS DATA ----------------------------------------------------------------------------------------------
  # Select ssb and log transform
  bsai.ssb <- bsai.ts %>%
    dplyr::select(YEAR, grep("ssb", colnames(.))) %>%
    mutate(bsai.ebs.pol.ssb = as.numeric(bsai.ebs.pol.ssb)) %>%
    pivot_longer(., !YEAR, names_to = "Stock", values_to="SSB") %>%
    mutate(SSB = log(SSB), Region = "BSAI") %>%
    rename(Year = YEAR)

ggplot(bsai.ssb, aes(x = Year, y = SSB))+
  geom_line(linewidth = 1)+
  facet_wrap(~Stock, scales = "free_y")
