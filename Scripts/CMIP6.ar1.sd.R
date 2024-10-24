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
    
    # Detrend data
    detrend.dat <- loess(annual.unsmoothed ~ year, TS.dat, span = 0.25, degree = 1)
    
    # Extract residuals
    TS.dat <- data.frame(year = TS.dat$year, annual.unsmoothed = detrend.dat$residuals)
    
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
    sd <-  rollapply(TS.dat$annual.unsmoothed, width = width, FUN = sd, fill = NA)
    
    # Make data frame of sd-cv
    data.frame(year = unique(TS.dat$year), sd = sd) -> win.dat
    
    # Calculate windows
    win.yr <- na.omit(win.dat) %>% pull(year)
    
    # Make data frame of ar1
    data.frame(year = win.yr, ar1 = ar1) -> ar1.dat
    
    # Join
    left_join(win.dat, ar1.dat) -> win.dat
    
    
    
    # Compile output
    sum.out.245 <- rbind(sum.out.245, cbind(win.dat, data.frame(model = models[ii],
                                                              region = regions[jj],
                                                              normalized_weight = norm_wt,
                                                              experiment = "hist_ssp245")))
    
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
    
    # Detrend data
    detrend.dat <- loess(annual.unsmoothed ~ year, TS.dat, span = 0.25, degree = 1)
    
    # Extract residuals
    TS.dat <- data.frame(year = TS.dat$year, annual.unsmoothed = detrend.dat$residuals)
    
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
    sd <-  rollapply(TS.dat$annual.unsmoothed, width = width, FUN = sd, fill = NA)
    
    # Make data frame of sd-cv
    data.frame(year = unique(TS.dat$year), sd = sd) -> win.dat
    
    # Calculate windows
    win.yr <- na.omit(win.dat) %>% pull(year)
    
    # Make data frame of ar1
    data.frame(year = win.yr, ar1 = ar1) -> ar1.dat
    
    # Join
    left_join(win.dat, ar1.dat) -> win.dat
    
    
   
    # Compile output
    sum.out.PI <- rbind(sum.out.PI, cbind(win.dat, data.frame(model = models[ii],
                                         region = regions[jj],
                                         normalized_weight = norm_wt,
                                         experiment = "piControl")))
    
  }
}

# BIND AND PLOT ---------------------
# Bind both ssp245 and pI control summaries
rbind(sum.out.PI, sum.out.245) -> sum.out

# Calculate weighted averages
sum.out %>%
  group_by(year, region, experiment) %>%
  reframe(ar1_wtdavg = weighted.mean(ar1, normalized_weight),
          sd_wtdavg = weighted.mean(sd, normalized_weight),
          ar1_wtdsd = sqrt(Hmisc::wtd.var(ar1, normalized_weight)),
          sd_wtdsd = sqrt(Hmisc::wtd.var(sd, normalized_weight))) -> wtd.sum.out

wtd.sum.out %>%
  filter(region == "Eastern_Bering_Sea", experiment == "hist_ssp245") %>%
  na.omit() -> dat.1

tau.ar1 = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.ar1.p = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$p.value, 2)
tau.sd = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.sd.p = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$p.value, 2)

df.1 <- data.frame(region = "Eastern Bering Sea", experiment = "hist_ssp245",
                   tau.ar1 = tau.ar1, tau.ar1.p = tau.ar1.p, tau.sd = tau.sd, tau.sd.p = tau.sd.p)


wtd.sum.out %>%
  filter(region == "Eastern_Bering_Sea", experiment == "piControl") %>%
  na.omit() -> dat.1

tau.ar1 = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.ar1.p = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$p.value, 2)
tau.sd = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.sd.p = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$p.value, 2)

df.2 <- data.frame(region = "Eastern Bering Sea", experiment = "piControl",
                   tau.ar1 = tau.ar1, tau.ar1.p = tau.ar1.p, tau.sd = tau.sd, tau.sd.p = tau.sd.p)


wtd.sum.out %>%
  filter(region == "Gulf_of_Alaska", experiment == "hist_ssp245") %>%
  na.omit() -> dat.1

tau.ar1 = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.ar1.p = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$p.value, 2)
tau.sd = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.sd.p = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$p.value, 2)

df.3 <- data.frame(region = "Gulf_of_Alaska", experiment = "hist_ssp245",
                   tau.ar1 = tau.ar1, tau.ar1.p = tau.ar1.p, tau.sd = tau.sd, tau.sd.p = tau.sd.p)


wtd.sum.out %>%
  filter(region == "Gulf_of_Alaska", experiment == "piControl") %>%
  na.omit() -> dat.1

tau.ar1 = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.ar1.p = round(cor.test(dat.1$ar1_wtdavg, dat.1$year, method = "kendall")$p.value, 2)
tau.sd = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$estimate, 2)
tau.sd.p = round(cor.test(dat.1$sd_wtdavg, dat.1$year, method = "kendall")$p.value, 2)

df.4 <- data.frame(region = "Gulf_of_Alaska", experiment = "piControl",
                   tau.ar1 = tau.ar1, tau.ar1.p = tau.ar1.p, tau.sd = tau.sd, tau.sd.p = tau.sd.p)

tau.dat <- rbind(df.1, df.2, df.3, df.4)

rownames(tau.dat) <- NULL

tau.dat %>%
  mutate(tau.ar1.sig = BH2(tau.ar1.p, alph = 0.05)$BHSig,
         tau.ar1.padj = p.adjust(tau.ar1.p, method = "fdr"),
         tau.sd.sig = BH2(tau.sd.p, alph = 0.05)$BHSig,
         tau.sd.padj = p.adjust(tau.sd.p, method = "fdr")) %>%
  right_join(tau.dat, .) -> tau.dat2

## Plot
labs <- c("Eastern Bering Sea", "Gulf of Alaska")
names(labs) <- c("Eastern_Bering_Sea", "Gulf_of_Alaska")

ggplot(wtd.sum.out, aes(year, ar1_wtdavg, color = region, linetype = experiment)) +
  geom_line(size = 1.25) +
  facet_wrap(~region, nrow = 2, scales = "free_y", labeller=  labeller(region = labs)) +
  theme_bw() +
  #geom_ribbon(wtd.sum.out, mapping = aes(ymin = ar1_wtdavg - ar1_wtdsd, ymax = ar1_wtdavg + ar1_wtdsd, fill = region), color = NA, alpha = 0.5)+
  scale_x_continuous(breaks = seq(min(wtd.sum.out$year), max(wtd.sum.out$year), by = 20),
                     labels = seq(min(wtd.sum.out$year), max(wtd.sum.out$year), by = 20))+
  scale_color_manual(values = c("#00AFBA", "#C28600"), guide = "none")+
  scale_fill_manual(values = c("#00AFBA", "#C28600"), name = "")+
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Anthropogenic radiative forcing\n(SSP245)", "Preindustrial\ncontrol"), name = "")+
  theme_bw()+
  ggtitle("Autocorrelation")+
  xlab("Year")+
  ylab("Weighted average")+
  theme(legend.direction= "horizontal",
        legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16))

ggsave("./Figures/CMIP6SST.AR1.png", width = 7, height = 7,units = "in")

ggplot(wtd.sum.out, aes(year, sd_wtdavg, color = region, linetype = experiment)) +
  geom_line(size = 1.25) +
  facet_wrap(~region, nrow = 2, , scales = "free_y", labeller=  labeller(region = labs)) +
  theme_bw() +
  #geom_ribbon(wtd.sum.out, mapping = aes(ymin = ar1_wtdavg - ar1_wtdsd, ymax = ar1_wtdavg + ar1_wtdsd, fill = region), color = NA, alpha = 0.5)+
  scale_x_continuous(breaks = seq(min(wtd.sum.out$year), max(wtd.sum.out$year), by = 20),
                     labels = seq(min(wtd.sum.out$year), max(wtd.sum.out$year), by = 20))+
  scale_color_manual(values = c("#00AFBA", "#C28600"), guide = "none")+
  scale_fill_manual(values = c("#00AFBA", "#C28600"), name = "")+
  scale_linetype_manual(values = c("solid", "twodash"), labels = c("Anthropogenic radiative forcing\n(SSP245)", "Preindustrial\ncontrol"), name = "")+
  theme_bw()+
  ggtitle("Standard deviation (°C)")+
  xlab("Year")+
  ylab("Weighted average")+
  theme(legend.direction= "horizontal",
        legend.position = "bottom",
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        strip.text = element_text(size = 16),
        legend.text = element_text(size = 14),
        title = element_text(size = 16))

ggsave("./Figures/CMIP6SST.SD.png", width = 7, height = 7, units = "in")
  
