#Calculating weighted mean AR1 and SD on 15-year rolling windows for SST anomalies from CMIP6 models

#PURPOSE: To load, process, and visualize GOA and EBS TS and climate timeseries

### LOAD PACKAGES/DATA -------------------------------------------------------------------------------------------------------

source("./Scripts/ts.processing.R")

ssp <- read.csv("./Data/CMIP6.anomaly.time.series.csv") # read in SST anomaly data
ssp245 <- read.csv("./Data/CMIP6.anomaly.time.series.ssp245.csv") # read in SST anomaly data

wts <- read.csv("./Data/normalized_CMIP6_weights.csv") # read in model weights


### Calculate AR1 and SD -------------

# Extract model names
models <- unique(ssp245$model)

# Specify region
regions = c("Eastern_Bering_Sea", "Gulf_of_Alaska")

# Run for loop 
sum.out.245 <- data.frame()

for(jj in 1:length(regions)){
  for(ii in 1:length(models)){
    
    # Filter anomaly data by model and region
    ssp245 %>%
      filter(model == models[ii], region == regions[jj]) %>%
      na.omit() -> TS.dat
    
    # Filter models weights by model and region
    wts %>%
      filter(model == model[ii], region == regions[jj]) %>%
      dplyr::select(normalized_weight) %>%
      pull() -> norm_wt
    
    # Specify sliding window width
    width = 15
    
    # Calculate rolling window AR1
    ar1 <- sapply(rollapply(TS.dat$annual.unsmoothed, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
    
    # Calculate rolling window SD
    sd <-  rollapply(TS.dat$annual.unsmoothed, width = width, FUN = sd)
    
    # Calculate windows
    window <- seq(min(TS.dat$year)+(width-1), max(TS.dat$year), by = 1)
    
    
    # Compile output
    sum.out.245 <- rbind(sum.out.245, data.frame(model = models[ii],
                                         ar1 = ar1,
                                         sd = sd,
                                         window = window,
                                         region = regions[jj],
                                         normalized_weight = norm_wt,
                                         experiment = "hist_ssp245"))
    
  }
}

# PI CONTROL ------------------------
# Extract model names
models <- unique(ssp$model)

# Specify region
regions = c("Eastern_Bering_Sea", "Gulf_of_Alaska")

# Run for loop 
sum.out.PI <- data.frame()

for(jj in 1:length(regions)){
  for(ii in 1:length(models)){
    
    # Filter anomaly data by model and region
    ssp %>%
      filter(model == models[ii], region == regions[jj], experiment == "piControl") %>%
      na.omit() -> TS.dat
    
    # Filter models weights by model and region
    wts %>%
      filter(model == model[ii], region == regions[jj]) %>%
      dplyr::select(normalized_weight) %>%
      pull() -> norm_wt
    
    # Specify sliding window width
    width = 15
    
    # Calculate rolling window AR1
    ar1 <- sapply(rollapply(TS.dat$annual.unsmoothed, width = width, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2) 
    
    # Calculate rolling window SD
    sd <-  rollapply(TS.dat$annual.unsmoothed, width = width, FUN = sd)
    
    # Calculate windows
    window <- seq(min(TS.dat$year)+(width-1), max(TS.dat$year), by = 1)
    
    
    # Compile output
    sum.out.PI <- rbind(sum.out.PI, data.frame(model = models[ii],
                                         ar1 = ar1,
                                         sd = sd,
                                         window = window,
                                         region = regions[jj],
                                         normalized_weight = norm_wt,
                                         experiment = "piControl"))
    
  }
}

# BIND AND PLOT ---------------------
# Bind both ssp245 and pI control summaries
rbind(sum.out.PI, sum.out.245) -> sum.out

# Calculate weighted averages
sum.out %>%
  group_by(window, region, experiment) %>%
  reframe(ar1_wtdavg = weighted.mean(ar1, normalized_weight),
          sd_wtdavg = weighted.mean(sd, normalized_weight),
          ar1_wtdsd = sqrt(Hmisc::wtd.var(ar1, normalized_weight)),
          sd_wtdsd = sqrt(Hmisc::wtd.var(sd, normalized_weight))) -> wtd.sum.out

ggplot(wtd.sum.out, aes(window, ar1_wtdavg, color = region, linetype = experiment)) +
  geom_line(size = 1.25) +
  facet_wrap(~region + experiment, nrow = 2) +
  theme_bw() +
  geom_ribbon(wtd.sum.out, mapping = aes(ymin = ar1_wtdavg - ar1_wtdsd, ymax = ar1_wtdavg + ar1_wtdsd, fill = region), color = NA, alpha = 0.5)+
  scale_x_continuous(breaks = seq(min(wtd.sum.out$window), max(wtd.sum.out$window), by = 15),
                     labels = seq(min(wtd.sum.out$window), max(wtd.sum.out$window), by = 15))+
  scale_color_manual(values = c("#6A6DB7", "#A34242"))+
  scale_fill_manual(values = c("#6A6DB7", "#A34242"))+
  theme_bw()+
  ylab("AR1 weighted AVG")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16))

ggsave("./Figures/CMIP6SST.AR1.png", width = 11, height = 8.5, units = "in")

ggplot(wtd.sum.out, aes(window, sd_wtdavg, color = region, linetype = experiment)) +
  geom_line(size = 1.25) +
  facet_wrap(~region, nrow = 2, scales = "free_y") +
  theme_bw() +
  scale_x_continuous(breaks = seq(min(wtd.sum.out$window), max(wtd.sum.out$window), by = 15),
                     labels = seq(min(wtd.sum.out$window), max(wtd.sum.out$window), by = 15))+
  scale_color_manual(values = c("#6A6DB7", "#A34242"))+
  theme_bw()+
  ylab("SD weighted AVG")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16))

ggsave("./Figures/CMIP6SST.SD.png", width = 11, height = 8.5, units = "in")
  
