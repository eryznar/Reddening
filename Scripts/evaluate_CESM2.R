# PURPOSE: to download CESM2 fcm and mdm SST/SLP model outputs and calculate ar1/SD

# Note: steps 1) and 2) take some time to process. If you just need to modify the final outputs/figures, proceed
# to step 3)

# Author: Emily Ryznar

# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# Get map layers
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

yrs <- 1850:2013 # these are years that are similar across SLP and SST models


# 1) DOWNLOAD CESM2 MODEL OUTPUTS ----------------------------------------------
# Authenticate google account
drive_auth()

# Identify drive folders to download from
fcm.sst <- drive_get(as_id("https://drive.google.com/drive/folders/1lycmXi992PgLAzBoSOWd3qommttHb6sc")) # FCM SST
mdm.sst <- drive_get(as_id("https://drive.google.com/drive/folders/1PbrbPKFzJIIVeCT6nGr0w9h-zoiCKnJQ")) # MDM SST

fcm.slp <- drive_get(as_id("https://drive.google.com/drive/folders/1cJh7VuEA1OhnVVgRAL_53y2c3XH7vIPZ")) # FCM SLP
mdm.slp <- drive_get(as_id("https://drive.google.com/drive/folders/1SU7zm-etGFzCgpNYfJ-jbW7jhKUgrF4K")) # MDM SLP

# Identify folders to download files into
fcm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/FCM/") #FCM SST
mdm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/MDM/") #MDM SST

fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP

# List the files in each folder
fcm.sst.files <- drive_ls(fcm.sst) # FCM SST
mdm.sst.files <- drive_ls(mdm.sst) # MDM SST

fcm.slp.files <- drive_ls(fcm.slp) # FCM SLP
mdm.slp.files <- drive_ls(mdm.slp) # MDM SLP

# Download the FCM SST files
files <- fcm.sst.files 
dir2 <- fcm.sst.dir

for (ii in 1:nrow(files)){
  print(paste0("Downloading file ", row(files)[ii],"/", nrow(files)))
  drive_download(as_id(files$id[ii]), path = paste0(dir2, files$name[ii]), overwrite = TRUE)
}

# Download the MDM SST files
files <- mdm.sst.files 
dir2 <- mdm.sst.dir

for (ii in 1:nrow(files)){
  print(paste0("Downloading file ", row(files)[ii],"/", nrow(files)))
  drive_download(as_id(files$id[ii]), path = paste0(dir2, files$name[ii]), overwrite = TRUE)
}

# Download the FCM SLP files
files <- fcm.slp.files 
dir2 <- fcm.slp.dir

for (ii in 1:nrow(files)){
  print(paste0("Downloading file ", row(files)[ii],"/", nrow(files)))
  drive_download(as_id(files$id[ii]), path = paste0(dir2, files$name[ii]), overwrite = TRUE)
}

# Download the MDM SLP files
files <- mdm.slp.files 
dir2 <- mdm.slp.dir

for (ii in 1:nrow(files)){
  print(paste0("Downloading file ", row(files)[ii],"/", nrow(files)))
  drive_download(as_id(files$id[ii]), path = paste0(dir2, files$name[ii]), overwrite = TRUE)
}

# 2) LOAD, STACK, and PROCESS CESM MODEL OUTPUTS --------------------------------
# Identify folders to download files from
fcm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/FCM/") #FCM SST
mdm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/MDM/") #MDM SST

fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP

# Extract time info for processing below (same across files)
files <- list.files(fcm.sst.dir, full.names = TRUE)

time_units <- ncmeta::nc_atts(files[1], "time") %>% 
  filter(name == "units") %>%
  pull(value)

unit_parts <- str_split(time_units, " since ")[[1]]
time_unit <- unit_parts[1]
origin_date <- ymd_hms(unit_parts[2])

  # FCM SST ----  
  files <- list.files(fcm.sst.dir, full.names = TRUE)
  fcm.sst <- tibble()
  
  files <- files[c(1:48, 50)]
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
                   #time = time > 711475) %>% # greater than 1947
      activate("SST") %>%
      hyper_tibble() %>%
      mutate(time = origin_date + lubridate::days(time),
             year = lubridate::year(time),
             month = lubridate::month(time),
             member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SSTa = SST - mean.month.SST,
             win.year = case_when((month %in% 11:12) ~ year + 1,
                                  TRUE ~ year)) %>% # compute anomalies
      group_by(lon, lat, year, member) %>% 
      reframe(mean.SSTa = mean(SSTa), 
              mean.winSSTa = mean(SSTa[month %in% c(11:12, 1:3)]))-> out # calculate mean annual SSTa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table

    out[, detrended_SSTa := { # output data table column
      fit <- lm(mean.SSTa ~ year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]
    
    
    # Calculate rolling window AR1 within data.table
    out[, ar1.SSTa := frollapply(detrended_SSTa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]
    
    # Calculate rolling window SD within data.table
    out[, sd.SSTa := frollapply(detrended_SSTa, n = 15, FUN = sd, 
                                    align = "center"), by = c("lon", "lat", "member")]
    
    # stack processed files
    fcm.sst <- bind_rows(fcm.sst, out)
    
  }
  
  setDT(fcm.sst)
  
  # Save file
  saveRDS(fcm.sst, paste0(dir, "Output/FCM_SSTa_ar1sd.rda"))
  
  # MDM SST ----
  files <- list.files(mdm.sst.dir, full.names = TRUE)
  mdm.sst <- tibble()
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
      #time = time > 711475) %>% # greater than 1947
      activate("SST") %>%
      hyper_tibble() %>%
      mutate(time = origin_date + lubridate::days(time),
             year = lubridate::year(time),
             month = lubridate::month(time),
             member = substr(files[ii], 68, 78),
             member= gsub(".postP.", ".", member)) %>% # extracting ensemble member #
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SSTa = SST - mean.month.SST,
             win.year = case_when((month %in% 11:12) ~ year + 1,
                                  TRUE ~ year)) %>% # compute anomalies
      group_by(lon, lat, year, member) %>% 
      reframe(mean.SSTa = mean(SSTa), 
              mean.winSSTa = mean(SSTa[month %in% c(11:12, 1:3)]))-> out # calculate mean annual SSTa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SSTa := { # output data table column
      fit <- lm(mean.SSTa ~ year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]
    
    
    # Calculate rolling window AR1 within data.table
    out[, ar1.SSTa := frollapply(detrended_SSTa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]
    
    # Calculate rolling window SD within data.table
    out[, sd.SSTa := frollapply(detrended_SSTa, n = 15, FUN = sd, 
                                align = "center"), by = c("lon", "lat", "member")]
    
    # stack processed files
    mdm.sst <- bind_rows(mdm.sst, out)
    
    
  }
  
  setDT(mdm.sst)
  
  # Save file
  saveRDS(mdm.sst, paste0(dir, "Output/MDM_SSTa_ar1sd.rda"))
  
  # FCM SLP ----
  files <- list.files(fcm.slp.dir, full.names = TRUE)
  fcm.slp <- tibble()
  
  files <- files[c(1:48, 50)]
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 191 & lon <= 209,
                   lat = lat >= 43 & lat <= 56) %>% # high activity AL area
      #time = time > 711475) %>% # greater than 1947
      activate("PSL") %>%
      hyper_tibble() %>%
      mutate(year = lubridate::year(time), # already in correct format
             month = lubridate::month(time),
             day = lubridate::day(time),
             member = substr(files[ii], 78, 85)) %>% # extracting ensemble member #
      filter(month %in% c(11:12, 1:3)) %>% #filter by winter months
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SLP = mean(PSL)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SLPa = PSL - mean.month.SLP, # compute anomalies
             win.year = case_when((month %in% 11:12) ~ year + 1, # calculate winter year
                                  TRUE ~ year)) %>% 
      filter(win.year %in% yrs) %>%
      group_by(lon, lat, win.year, member) %>% 
      reframe(mean.SLPa = mean(SLPa)) -> out # calculate mean annual SLPa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SLPa := { # output column name
      fit <- lm(mean.SLPa ~ win.year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]

    # Calculate rolling window AR1 within data.table
    out[, ar1.SLPa := frollapply(detrended_SLPa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]
    
    # Calculate rolling window SD within data.table
    out[, sd.SLPa := frollapply(detrended_SLPa, n = 15, FUN = sd, 
                                align = "center"), by = c("lon", "lat", "member")]
    
    # stack processed files
    fcm.slp <- bind_rows(fcm.slp, out)
    
    
  }
  
  setDT(fcm.slp)
  
  # Save file
  saveRDS(fcm.slp, paste0(dir, "Output/FCM_winterSLPa_ar1sd.rda"))
  
  # MDM SLP ----
  files <- list.files(mdm.slp.dir, full.names = TRUE)
  mdm.slp <- tibble()
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 191 & lon <= 209,
                   lat = lat >= 43 & lat <= 56) %>% # high activity AL area
      #time = time > 711475) %>% # greater than 1947
      activate("PSL") %>%
      hyper_tibble() %>%
      mutate(year = lubridate::year(time), # already in correct format
             month = lubridate::month(time),
             day = lubridate::day(time),
             member = substr(files[ii], 68, 72)) %>% # extracting ensemble member #
      filter(month %in% c(11:12, 1:3),
             year %in% yrs) %>% #filter by winter months
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SLP = mean(PSL)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SLPa = PSL - mean.month.SLP, # compute anomalies
             win.year = case_when((month %in% 11:12) ~ year + 1, # calculate winter year
                                  TRUE ~ year)) %>%
      filter(win.year %in% yrs) %>%
      group_by(lon, lat, win.year, member) %>% 
      reframe(mean.SLPa = mean(SLPa)) -> out # calculate mean annual SLPa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SLPa := { # output column name
      fit <- lm(mean.SLPa ~ win.year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]
    
    # Calculate rolling window AR1 within data.table
    out[, ar1.SLPa := frollapply(detrended_SLPa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]
    
    # Calculate rolling window SD within data.table
    out[, sd.SLPa := frollapply(detrended_SLPa, n = 15, FUN = sd, 
                                align = "center"), by = c("lon", "lat", "member")]
    
    # stack processed files
    mdm.slp <- bind_rows(mdm.slp, out)
    
    
  }
  
  setDT(mdm.slp)
  
  # Save file
  saveRDS(mdm.slp, paste0(dir, "Output/MDM_winterSLPa_ar1sd.rda"))
  
  
# 3) CALCULATE/PLOT CELL-WISE SST AR1 SD, SST SD, and MEAN AR1 ACROSS ENSEMBLE ------------------------
  # FCM SST ----
  fcm.sst <- readRDS(paste0(dir, "Output/FCM_SSTa_ar1sd.rda"))

  # Calculate grid cell AR1
  fcm.sst %>%
    na.omit() %>%
    group_by(lon, lat, member) %>%
    reframe(ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
            ar1.cv = sd(ar1.SSTa)/(abs(mean(ar1.SSTa))), # sd of AR1 of rolling window timeseries
            ar1.mean = mean(ar1.SSTa), # mean AR1 of rolling window timeseries
            sd.sd = sd(sd.SSTa), # sd of SD of rolling window timeseries
            sd.mean = mean(sd.SSTa)) %>% # mean SD of rolling window timeseries
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat)) %>%
    distinct() -> fcm.ar1.sd.dat

  # MDM SST ----
  mdm.sst <- readRDS(paste0(dir, "Output/MDM_SSTa_ar1sd.rda"))

  # Calculate grid cell AR1
  mdm.sst %>%
    na.omit() %>%
    group_by(lon, lat, member) %>%
    reframe(ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
            ar1.cv = sd(ar1.SSTa)/(abs(mean(ar1.SSTa))), # sd of AR1 of rolling window timeseries
            ar1.mean = mean(ar1.SSTa), # mean AR1 of rolling window timeseries
            sd.sd = sd(sd.SSTa), # sd of SD of rolling window timeseries
            sd.mean = mean(sd.SSTa)) %>% # mean SD of rolling window timeseries
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat)) %>%
    distinct() -> mdm.ar1.sd.dat



# Calculate mean SD AR1 and AR1 by grid cell across ensemble members
fcm.ar1.sd.dat %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.cv = mean(ar1.cv),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Fully coupled") -> fcm.ar1.sd.dat2


mdm.ar1.sd.dat %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.cv = mean(ar1.cv),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Mechanically decoupled") -> mdm.ar1.sd.dat2

# Join
plot.dat <- rbind(fcm.ar1.sd.dat2, mdm.ar1.sd.dat2)

# Calculate differences
plot.dat %>%
  pivot_wider(names_from = type, values_from = c(ar1.sd, ar1.cv, ar1.mean, sd.mean)) %>%
  mutate(ar1.sd.diff = scale(`ar1.sd_Fully coupled` - `ar1.sd_Mechanically decoupled`)[,1],
         ar1.cv.diff = scale(`ar1.cv_Fully coupled` - `ar1.cv_Mechanically decoupled`)[,1],
         ar1.mean.diff = scale(`ar1.mean_Fully coupled` - `ar1.mean_Mechanically decoupled`)[,1],
         sd.mean.diff = scale(`sd.mean_Fully coupled` - `sd.mean_Mechanically decoupled`)[,1],
         type = "Difference") %>%
  dplyr::select(lon, lat, ar1.sd.diff, ar1.cv.diff, ar1.mean.diff, sd.mean.diff, type) -> diff.dat


# Plot AR1 SD and difference
ggplot()+
  geom_tile(plot.dat, mapping= aes(lon, lat, fill = ar1.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type, nrow = 2)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  median(plot.dat$ar1.sd),
                       name = "AR1 sd")+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> ar1.sd.plot

ggplot()+
  geom_tile(diff.dat, mapping= aes(lon, lat, fill = ar1.sd.diff))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  0,
                       name = "FCM-MDM diff",
                       limits = c(max(abs(diff.dat$ar1.sd.diff))*-1, max(abs(diff.dat$ar1.sd.diff))))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> ar1.sd.diffplot


ar1.sd.plot + ar1.sd.diffplot + plot_layout(nrow = 2, ncol = 1, byrow = TRUE,
                                            widths = c(1, 1), heights = c(1, 0.5),
                                            axes = "collect")

ggsave("./Figures/CESM2_AR1_SD.png", height= 7, width = 5, units = "in")

# Plot AR1 mean and difference
ggplot()+
  geom_tile(plot.dat, mapping= aes(lon, lat, fill = ar1.mean))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type, nrow = 2)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  median(plot.dat$ar1.mean),
                       name = "mean AR1")+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> ar1.mean.plot

ggplot()+
  geom_tile(diff.dat, mapping= aes(lon, lat, fill = ar1.mean.diff))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  0,
                       name = "FCM-MDM diff",
                       limits = c(max(abs(diff.dat$ar1.mean.diff))*-1, max(abs(diff.dat$ar1.mean.diff))))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> ar1.mean.diffplot


ar1.mean.plot + ar1.mean.diffplot + plot_layout(nrow = 2, ncol = 1, byrow = TRUE,
                                            widths = c(1, 1), heights = c(1, 0.5),
                                            axes = "collect")

ggsave("./Figures/CESM2_AR1_MEAN.png", height= 7, width = 5, units = "in")

# Plot SD mean and difference
ggplot()+
  geom_tile(plot.dat, mapping= aes(lon, lat, fill = sd.mean))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type, nrow = 2)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  median(plot.dat$sd.mean),
                       name = "mean SD")+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> sd.mean.plot

ggplot()+
  geom_tile(diff.dat, mapping= aes(lon, lat, fill = sd.mean.diff))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                       midpoint=  0,
                       name = "FCM-MDM diff",
                       limits = c(max(abs(diff.dat$sd.mean.diff))*-1, max(abs(diff.dat$sd.mean.diff))))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> sd.mean.diffplot


sd.mean.plot + sd.mean.diffplot + plot_layout(nrow = 2, ncol = 1, byrow = TRUE,
                                                widths = c(1, 1), heights = c(1, 0.5),
                                                axes = "collect")

ggsave("./Figures/CESM2_SD_MEAN.png", height= 7, width = 5, units = "in")


# # 4) EVALUATE PERIODS OF HIGH SLP VARIABILITY ----
# # slp fcm
# slp.fcm <- readRDS(paste0(dir, "Output/FCM_winterSLPa_ar1sd.rda"))
# 
# mems <- unique(slp.fcm$member)
# 
# 
# slp.fcm %>%
#   #filter(member == mems[16]) %>% # isolating one ensemble member
#   group_by(win.year) %>%
#   reframe(sd = mean(sd.SLPa)) %>%
#   na.omit()-> tt
# 
# # Seems like 1984-1921 has high AL variability and 1922-1947 has low
# ggplot(tt %>% filter(win.year > 1857), aes(win.year, sd))+
#   geom_point()+
#   geom_line() + 
#   scale_x_continuous(breaks = seq(min(tt$win.year), max(tt$win.year), by = 10))+
#   theme_bw() +
#   ylab("SLP standard deviation")+
#   ggtitle("AL high activity area winter SLPa")
# 
# #slp mdm
# slp.mdm <- readRDS(paste0(dir, "Output/mdm_winterSLPa_ar1sd.rda"))
# 
# mems <- unique(slp.mdm$member)
# 
# 
# slp.mdm %>%
#   #filter(member == mems[16]) %>% # isolating one ensemble member
#   group_by(win.year) %>%
#   reframe(sd = mean(sd.SLPa)) %>%
#   na.omit()-> tt
# 
# # Seems like 1984-1921 has high AL variability and 1922-1947 has low
# ggplot(tt %>% filter(win.year > 1857), aes(win.year, sd))+
#   geom_point()+
#   geom_line() + 
#   scale_x_continuous(breaks = seq(min(tt$win.year), max(tt$win.year), by = 10))+
#   theme_bw() +
#   ylab("SLP standard deviation")+
#   ggtitle("AL high activity area winter SLPa")
# 
# #1916-1935 low
# #1966-1990 high
# 
# length(1916:1935) #low
# length(1966:1985) #high
# # # By ensemble member
# # slp.fcm %>%
# #   group_by(win.year, member) %>%
# #   reframe(sd = mean(sd.SLPa)) %>%
# #   na.omit()-> tt
# # 
# # 
# # ggplot(tt %>% filter(win.year > 1857, member %in% mems[1:10]), aes(win.year, sd, group = member))+
# #   #geom_point()+
# #   geom_line() + 
# #   scale_x_continuous(breaks = seq(min(tt$win.year), max(tt$win.year), by = 10))+
# #   theme_bw() +
# #   ylab("SLP standard deviation")+
# #   ggtitle("AL high activity area winter SLPa")+
# #   theme(legend.position = "none")
# 
# highyrs <- 1966:1985
# lowyrs <- 1916:1935
# 
#   # Process FCM and MDM SST in winter for comparison in high and low SLP var periods ----
#   files <- list.files(fcm.sst.dir, full.names = TRUE)
#   
#   time_units <- ncmeta::nc_atts(files[1], "time") %>% 
#     filter(name == "units") %>%
#     pull(value)
#   
#   unit_parts <- str_split(time_units, " since ")[[1]]
#   time_unit <- unit_parts[1]
#   origin_date <- ymd_hms(unit_parts[2])
#   
#   # FCM SST  
#   files <- list.files(fcm.sst.dir, full.names = TRUE)
#   
#   fcm.win.sst2 <- tibble()
#   
#   for(ii in 1:length(files)){
#     
#     print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
#     
#     # load and process for high var periods
#     tidync(files[ii]) %>%
#       hyper_filter(lon = lon >= 125 & lon <= 255,
#                    lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
#       #time = time > 711475) %>% # greater than 1947
#       activate("SST") %>%
#       hyper_tibble() %>%
#       mutate(time = origin_date + lubridate::days(time),
#              year = lubridate::year(time),
#              month = lubridate::month(time),
#              member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
#       group_by(lat, lon, month, member) %>%
#       mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
#       ungroup() %>%
#       mutate(SSTa = SST - mean.month.SST,
#              win.year = case_when((month %in% 11:12) ~ year + 1,
#                                   TRUE ~ year)) %>% # compute anomalies
#       filter(win.year %in% c(highyrs, lowyrs),
#              month %in% c(11:12, 1:3)) %>%
#       mutate(period = case_when((win.year %in% highyrs) ~ "high",
#                                 TRUE ~ "low")) %>%
#       group_by(lon, lat, win.year, member, period) %>% 
#       reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell
#     
#     
#     # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
#     setDT(out) # convert to data.table
#     
#     out[, detrended_winSSTa := { # output data table column
#       fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
#       residuals(fit)           # Extract residuals as detrended values
#     }, by = .(lat, lon, period)]
#     
#     
#     # Calculate rolling window AR1 within data.table
#     out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
#       acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
#       return(acf_result$acf[2])
#     }, align = "center"), by = c("lon", "lat", "member", "period")]
#     
#     # Calculate rolling window SD within data.table
#     out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd, 
#                                 align = "center"), by = c("lon", "lat", "member", "period")]
#     
#     # stack processed files
#     fcm.win.sst2 <- bind_rows(fcm.win.sst2, out)
#     
#   }
#   
#   setDT(fcm.win.sst)
#   
#   # Save file
#   saveRDS(fcm.win.sst, paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda"))
#   
#   # MDM SST  
#   files <- list.files(mdm.sst.dir, full.names = TRUE)
#   
#   mdm.win.sst <- tibble()
#   
#   for(ii in 1:length(files)){
#     
#     print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
#     
#     # load and process for high var periods
#     tidync(files[ii]) %>%
#       hyper_filter(lon = lon >= 125 & lon <= 255,
#                    lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
#       #time = time > 711475) %>% # greater than 1947
#       activate("SST") %>%
#       hyper_tibble() %>%
#       mutate(time = origin_date + lubridate::days(time),
#              year = lubridate::year(time),
#              month = lubridate::month(time),
#              member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
#       group_by(lat, lon, month, member) %>%
#       mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
#       ungroup() %>%
#       mutate(SSTa = SST - mean.month.SST,
#              win.year = case_when((month %in% 11:12) ~ year + 1,
#                                   TRUE ~ year)) %>% # compute anomalies
#       filter(win.year %in% c(highyrs, lowyrs),
#              month %in% c(11:12, 1:3)) %>%
#       mutate(period = case_when((win.year %in% highyrs) ~ "high",
#                                 TRUE ~ "low")) %>%
#       group_by(lon, lat, win.year, member, period) %>% 
#       reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell
#     
#     
#     # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
#     setDT(out) # convert to data.table
#     
#     out[, detrended_winSSTa := { # output data table column
#       fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
#       residuals(fit)           # Extract residuals as detrended values
#     }, by = .(lat, lon, period)]
#     
#     
#     # Calculate rolling window AR1 within data.table
#     out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
#       acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
#       return(acf_result$acf[2])
#     }, align = "center"), by = c("lon", "lat", "member", "period")]
#     
#     # Calculate rolling window SD within data.table
#     out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd, 
#                                    align = "center"), by = c("lon", "lat", "member", "period")]
#     
#     # stack processed files
#     mdm.win.sst <- bind_rows(mdm.win.sst, out)
#     
#   }
#   
#   setDT(mdm.win.sst)
#   
#   # Save file
#   saveRDS(mdm.win.sst, paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda"))
#   
#   
#   # Calculate and plot cell-wise AR1, SD, mean AR1 ----
#   fcm.sst <- readRDS(paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda"))
#   
#   # Calculate grid cell AR1
#   fcm.sst %>%
#     na.omit() %>%
#     group_by(lon, lat, member, period) %>%
#     reframe(ar1.sd = sd(ar1.winSSTa), # sd of AR1 of rolling window timeseries
#             ar1.cv = (sd(ar1.winSSTa)/mean(ar1.winSSTa))*100, # sd of AR1 of rolling window timeseries
#             ar1.mean = mean(ar1.winSSTa), # mean AR1 of rolling window timeseries
#             sd.mean = mean(sd.winSSTa)) %>% # mean SD of rolling window timeseries
#     mutate(lon = as.numeric(lon),
#            lat = as.numeric(lat)) %>%
#     distinct() -> fcm.ar1.sd.dat
#   
#   
#   # Calculate grid cell AR1
#   mdm.sst %>%
#     na.omit() %>%
#     group_by(lon, lat, member, period) %>%
#     reframe(ar1.sd = sd(ar1.winSSTa), # sd of AR1 of rolling window timeseries
#             ar1.cv = (sd(ar1.winSSTa)/mean(ar1.winSSTa))*100, # sd of AR1 of rolling window timeseries
#             ar1.mean = mean(ar1.winSSTa), # mean AR1 of rolling window timeseries
#             sd.mean = mean(sd.winSSTa)) %>% # mean SD of rolling window timeseries
#     mutate(lon = as.numeric(lon),
#            lat = as.numeric(lat)) %>%
#     distinct() -> mdm.ar1.sd.dat
#   
#   
#   
#   # Calculate mean SD AR1 and AR1 by grid cell across ensemble members
#   fcm.ar1.sd.dat %>%
#     group_by(lon, lat, period) %>%
#     reframe(ar1.sd = mean(ar1.sd),
#             ar1.cv = mean(ar1.cv),
#             ar1.mean = mean(ar1.mean),
#             sd.mean = mean(sd.mean)) %>%
#     mutate(type = "FCM") -> fcm.ar1.sd.dat2
#   
#   
#   mdm.ar1.sd.dat %>%
#     group_by(lon, lat, period) %>%
#     reframe(ar1.sd = mean(ar1.sd),
#             ar1.cv = mean(ar1.cv),
#             ar1.mean = mean(ar1.mean),
#             sd.mean = mean(sd.mean)) %>%
#     mutate(type = "MDM") -> mdm.ar1.sd.dat2
#   
#   # Join
#   plot.dat <- rbind(fcm.ar1.sd.dat2, mdm.ar1.sd.dat2) %>%
#     mutate(period = case_when((period == "high") ~ "High",
#                               TRUE ~ "Low"))
#   
#   # Pivot longer and calculate differences
#   diff.dat <- pivot_longer(cols = c(4:7), values_to = "value", names_to = "name", data = plot.dat) %>%
#               group_by(lon, lat, name, period) %>%
#               mutate(fcm.mdm_diff = value[type == "FCM"] - value[type == "MDM"]) %>%
#               ungroup() %>%
#               group_by(name, period) %>%
#               mutate(scale_diff = scale(fcm.mdm_diff)[,1]) %>% 
#               ungroup() %>%
#               mutate(type2 = "Difference")
#   
#   # Plot AR1 SD and difference
#   dat <- diff.dat %>% filter(name == "ar1.sd")
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = value))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type~period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  median(dat$value),
#                          name = "AR1 sd")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           axis.text.x = element_blank()) -> ar1.sd.plot
#   
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = scale_diff))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type2 ~ period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  0,
#                          name = "FCM-MDM diff",
#                          limits = c(max(abs(dat$scale_diff))*-1, max(abs(dat$scale_diff))))+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           strip.text.x = element_blank()) -> ar1.sd.diffplot
#   
#   
#   ar1.sd.plot + ar1.sd.diffplot + plot_layout(nrow = 2, ncol = 1, byrow = TRUE, 
#                                               widths = c(1, 1), heights = c(1, 0.5),
#                                               axes = "collect")
#   
#   ggsave("./Figures/CESM2_AR1_SD_ALhighlow.png", height= 7, width = 7, units = "in")
#   
#   # Plot AR1 mean and difference
#   dat <- diff.dat %>% filter(name == "ar1.mean")
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = value))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type~period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  median(dat$value),
#                          name = "AR1 mean")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           axis.text.x = element_blank()) -> plot.1
#   
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = scale_diff))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type2 ~ period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  0,
#                          name = "FCM-MDM diff",
#                          limits = c(max(abs(dat$scale_diff))*-1, max(abs(dat$scale_diff))))+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           strip.text.x = element_blank()) -> plot.2
#   
#   
#   plot.1 + plot.2 + plot_layout(nrow = 2, ncol = 1, byrow = TRUE, 
#                                               widths = c(1, 1), heights = c(1, 0.5),
#                                               axes = "collect")
#   
#   ggsave("./Figures/CESM2_AR1_MEAN_ALhighlow.png", height= 7, width = 7, units = "in")
#   
#   
#   # Plot SD mean and difference
#   dat <- diff.dat %>% filter(name == "sd.mean")
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = value))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type~period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  median(dat$value),
#                          name = "SD mean")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           axis.text.x = element_blank()) -> plot.1
#   
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = scale_diff))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type2 ~ period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  0,
#                          name = "FCM-MDM diff",
#                          limits = c(max(abs(dat$scale_diff))*-1, max(abs(dat$scale_diff))))+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           strip.text.x = element_blank()) -> plot.2
#   
#   
#   plot.1 + plot.2 + plot_layout(nrow = 2, ncol = 1, byrow = TRUE, 
#                                               widths = c(1, 1), heights = c(1, 0.5),
#                                               axes = "collect")
#   
#   ggsave("./Figures/CESM2_SD_MEAN_ALhighlow.png", height= 7, width = 7, units = "in")
#   
#   # Plot AR1 CV and difference
#   dat <- diff.dat %>% filter(name == "ar1.cv")
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = value))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type~period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  median(dat$value),
#                          name = "AR1 cv")+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           axis.text.x = element_blank()) -> plot.1
#   
#   
#   ggplot()+
#     geom_tile(dat, mapping= aes(lon, lat, fill = scale_diff))+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
#     facet_grid(type2 ~ period)+
#     xlab("Latitude")+
#     ylab("Longitude")+
#     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
#                          midpoint=  0,
#                          name = "FCM-MDM diff",
#                          limits = c(max(abs(dat$scale_diff))*-1, max(abs(dat$scale_diff))))+
#     theme_bw()+
#     theme(plot.title = element_text(size = 10),
#           legend.title = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           strip.text.x = element_blank()) -> plot.2
#   
#   
#   plot.1 + plot.2 + plot_layout(nrow = 2, ncol = 1, byrow = TRUE, 
#                                               widths = c(1, 1), heights = c(1, 0.5),
#                                               axes = "collect")
#   
#   ggsave("./Figures/CESM2_AR1_CV_ALhighlow.png", height= 7, width = 7, units = "in")
#   
#   
# 5) AL SLPa and WINTER SSTa REGRESSIONS ----
  # Process CESM2 SST for winter and all years ----
  files <- list.files(fcm.sst.dir, full.names = TRUE)

  time_units <- ncmeta::nc_atts(files[1], "time") %>%
    filter(name == "units") %>%
    pull(value)

  unit_parts <- str_split(time_units, " since ")[[1]]
  time_unit <- unit_parts[1]
  origin_date <- ymd_hms(unit_parts[2])

  # FCM SST
  files <- list.files(fcm.sst.dir, full.names = TRUE)
  files <- files[c(1:48, 50)] # 49 is throwing an error

  fcm.win.sst <- tibble()

  for(ii in 1:length(files)){

    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking

    # load and process for high var periods
    tidync(files[1]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
      #time = time > 711475) %>% # greater than 1947
      activate("SST") %>%
      hyper_tibble() %>%
      mutate(time = origin_date + lubridate::days(time),
             year = lubridate::year(time),
             month = lubridate::month(time),
             member = substr(files[1], 81, 88)) %>% # extracting ensemble member #
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SSTa = SST - mean.month.SST,
             win.year = case_when((month %in% 11:12) ~ year + 1,
                                  TRUE ~ year)) %>% # compute anomalies
      filter(month %in% c(11:12, 1:3),
             year %in% yrs) %>%
      group_by(lon, lat, win.year, member) %>%
      reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell


    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table

    out[, detrended_winSSTa := { # output data table column
      fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]


    # Calculate rolling window AR1 within data.table
    out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]

    # Calculate rolling window SD within data.table
    out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd,
                                   align = "center"), by = c("lon", "lat", "member")]

    # stack processed files
    fcm.win.sst <- bind_rows(fcm.win.sst, out)

  }

  setDT(fcm.win.sst)

  # Save file
  saveRDS(fcm.win.sst, paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda"))

  # MDM SST
  files <- list.files(mdm.sst.dir, full.names = TRUE)

  mdm.win.sst <- tibble()

  for(ii in 1:length(files)){

    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking


    # load and process for high var periods
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
      #time = time > 711475) %>% # greater than 1947
      activate("SST") %>%
      hyper_tibble() %>%
      mutate(time = origin_date + lubridate::days(time),
             year = lubridate::year(time),
             month = lubridate::month(time),
             member = substr(files[ii], 68, 78),
             member= gsub(".postP.", ".", member)) %>% # extracting ensemble member #
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SSTa = SST - mean.month.SST,
             win.year = case_when((month %in% 11:12) ~ year + 1,
                                  TRUE ~ year)) %>% # compute anomalies
      filter(month %in% c(11:12, 1:3),
             year %in% yrs) %>%
      group_by(lon, lat, win.year, member) %>%
      reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell


    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table

    out[, detrended_winSSTa := { # output data table column
      fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]


    # Calculate rolling window AR1 within data.table
    out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
      acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
      return(acf_result$acf[2])
    }, align = "center"), by = c("lon", "lat", "member")]

    # Calculate rolling window SD within data.table
    out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd,
                                   align = "center"), by = c("lon", "lat", "member")]

    # stack processed files
    mdm.win.sst <- bind_rows(mdm.win.sst, out)

  }

  setDT(mdm.win.sst)

  # Save file
  saveRDS(mdm.win.sst, paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda"))

  # Regressions using CESM output -----
    ## Read in CESM2 output files
    sst <- rbind(readRDS(paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda")) %>%
                 mutate(model = "FCM"),
                 readRDS(paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda")) %>%
                 mutate(model = "MDM")) %>%
           mutate(lat = as.numeric(lat), lon = as.numeric(lon))

    slp <- rbind(readRDS(paste0(dir, "Output/FCM_winterSLPa_ar1sd.rda")) %>%
                 mutate(model = "FCM"),
                 readRDS(paste0(dir, "Output/MDM_winterSLPa_ar1sd.rda")) %>%
                 mutate(model = "MDM")) %>%
           mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
          filter(member %in% sst$member)


    ## Filter sst to EBS and GOA
    # Specify EBS polygon
    ebs.x <- c(183, 183, 203, 203, 191)
    #ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
    ebs.y <- c(53, 65, 65, 57.5, 53)

    # Filter sst initially to make things go quicker
    sst.check <- sst %>%
                filter(lon >= 170 & lon <= 210,
                       lat >= 50 & lat <= 70)

    # Create EBS polygon, filter sst not within that polygon
    xp <- cbind(ebs.x, ebs.y)
    loc= sst.check %>% dplyr::select(lon, lat)
    check <- in.poly(loc, xp=xp)

    sst.ebs <- cbind(sst.check, check) %>%
                filter(check == TRUE) %>%
                dplyr::select(!check)

    # Plot to check
    sst.ebs %>%
      filter(win.year == 1850, member == "1301.020", model == "FCM") -> tt

    ggplot()+
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
      coord_cartesian(ylim = c(50, 67), xlim = c(180, 205), expand = FALSE)+
      geom_point(tt, mapping= aes(lon, lat))

    sst.ebs <- sst.ebs %>%
                na.omit() %>%
                group_by(win.year, member, model) %>%
                reframe(ar1.winSSTa = mean(ar1.winSSTa))

    # Specify GOA polygon
    goa.x <- c(201, 201, 205, 208, 225, 231, 201)
    #goa.x <- ifelse(goa.x > 180, goa.x-360, goa.x)
    goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

    # Filter sst initially to make things go quicker
    sst.check <- sst %>%
      filter(lon >= 190 & lon <= 250,
             lat >= 50 & lat <= 70)

    # Create EBS polygon, filter sst not within that polygon
    xp <- cbind(goa.x, goa.y)
    loc= sst.check %>% dplyr::select(lon, lat)
    check <- in.poly(loc, xp=xp)

    sst.goa <- cbind(sst.check, check) %>%
      filter(check == TRUE) %>%
      dplyr::select(!check)

    # Plot to check
    sst.goa %>%
      filter(win.year == 1850, member == "1301.020", model == "FCM") -> tt

    ggplot()+
      geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
      coord_cartesian(ylim = c(50, 67), xlim = c(190, 250), expand = FALSE)+
      geom_point(tt, mapping= aes(lon, lat))

    sst.goa <- sst.goa %>%
                na.omit() %>%
                group_by(win.year, member, model) %>%
                reframe(ar1.winSSTa = mean(ar1.winSSTa))

    ## Run regressions
    slp <- slp %>%
            na.omit() %>%
            group_by(win.year, member, model) %>%
            reframe(sd.winSLPa = mean(sd.SLPa))


    # EBS
    right_join(slp, sst.ebs) %>%
      na.omit()-> mod.dat

    mod.ebs <- lme(ar1.winSSTa ~ sd.winSLPa+model, data = mod.dat, random = ~ 1| member,
                   correlation = corAR1(form = ~ win.year|member))

    summary(mod.ebs)

    mod.ebs2 <- gamm(
      ar1.winSSTa ~ s(sd.winSLPa, bs = "cr") + model,  # Fixed effects
      random = list(member = ~1),      # Random effects (named list)
      correlation = corAR1(form = ~ win.year | member),  # AR(1) structure
      data = mod.dat
    )

    summary(mod.ebs2$gam)

    AICc(mod.ebs, mod.ebs2)

    # Effect plots
    rr <- Effect(c("sd.winSLPa", "model"), mod.ebs, partial.residuals = TRUE)
    plot(rr,  main = "EBS partial residuals", show_residuals = TRUE, show_residuals_line = FALSE)

    rr <- Effect(c("sd.winSLPa", "model"), mod.ebs, confidence.level = 0.95)
    plot(rr, main = "EBS fitted values with 95% CI")

    tt <- predictorEffect("sd.winSLPa", mod.ebs)

    # GOA
    right_join(slp, sst.goa) %>%
      na.omit()-> mod.dat

    mod.goa <- lme(ar1.winSSTa ~ sd.winSLPa+model, data = mod.dat, random = ~ 1| member,
                  correlation = corAR1(form = ~ win.year|member))

    summary(mod.goa)

    mod.goa2 <- gamm(
      ar1.winSSTa ~ s(sd.winSLPa)+model,  # Fixed effects
      random = list(member = ~1),      # Random effects (named list)
      correlation = corAR1(form = ~ win.year | member),  # AR(1) structure
      data = mod.dat
    )

    summary(mod.goa2$gam)

    AICc(mod.goa, mod.goa2)

    # Effect plots
    rr <- Effect(c("sd.winSLPa", "model"), mod.goa, partial.residuals = TRUE)
    plot(rr, partial.residuals = TRUE, main = "GOA partial residuals")

    rr <- Effect(c("sd.winSLPa", "model"), mod.goa, confidence.level = 0.95)
    plot(rr, main = "GOA fitted values with 95% CI")

  # Regressions using observation data ----
  ## Read in winter SST and SLP SD and AR1 data (data is detrended)
  model.dat <- read.csv(paste0(dir, "Output/winSSTSLP_regressiondat.csv"))

  ## Fit models
    # Model: EBS sst AR1 x slp AR1
    mod.1 <- gls(ebs.sst.ar1~slp.sd,
                 data= model.dat, correlation = corAR1(form = ~year))

    summary(mod.1)

    mod.2 <- gamm(ebs.sst.ar1~s(slp.sd, bs = "cr"), data = model.dat, correlation = corAR1(form = ~year))

    summary(mod.2$gam)

    AICc(mod.1, mod.2)

    # Model: GOA sst AR1 x slp AR1
    mod.1 <- gls(goa.sst.ar1~slp.sd,
                 data= model.dat, correlation = corAR1(form = ~year))

    summary(mod.1)

    mod.2 <-  gamm(goa.sst.ar1~s(slp.sd, bs = "cr"), data = model.dat, correlation = corAR1(form = ~year))

    summary(mod.2$gam)

    AICc(mod.1, mod.2)


# # 6) AL SLPa and WINTER SSTa REGRESSION RANDOMIZATIONS ----
#   # Observations ----
#   ## Read in winter SST and SLP SD and AR1 data (data is detrended)
#   model.dat <- read.csv(paste0(dir, "Output/winSSTSLP_regressiondat.csv"))
#     
#   ## Get AR1 and SD of SLPa SD
#   SLPaSD.sd <- sd(model.dat$slp.sd)
#   
#   SLPaSD.ar1 <- acf(model.dat$slp.sd, plot = FALSE, na.action = na.pass)$acf[2]
#  
#   ## Create simulation dataframe
#   sim.dat <- data.frame()
#   
#   ## Create randomized SLPa SD timeseries, fit models for GOA and EBS SSTa
#   iter <- 1:10000
#   
#   for(ii in 1:length(iter)){
#     # Simulate SLPa SD timeseries
#     SLPts <- arima.sim(model = list(order = c(1,0,0), 
#                            ar = SLPaSD.ar1, 
#                            sd= SLPaSD.sd),
#                            n = nrow(model.dat)) 
#     
#     # EBS
#     mod.ebs <- lm(model.dat$ebs.sst.ar1 ~ SLPts)
#     
#     # GOA
#     mod.goa <- lm(model.dat$goa.sst.ar1 ~ SLPts)
#     
#     # Create output df
#     out <- data.frame(iteration = iter[ii], 
#                       region = c("EBS", "GOA"), 
#                       F_stat = c(round(glance(mod.ebs)$statistic,2), round(glance(mod.goa)$statistic,2)),
#                       p_val = c(round(glance(mod.ebs)$p.value,2), round(glance(mod.goa)$p.value,2)))
#     
#     # Bind across iterations
#     sim.dat <- rbind(sim.dat, out)
#     
#   }
#   
#   
#   ## Fit models with observed SLP.sd, extract Fstat and pval
#   mod.ebs <- lm(ebs.sst.ar1~slp.sd, data= model.dat)
#   ebs.F <- glance(mod.ebs)$statistic
#   ebs.p <- glance(mod.ebs)$p.value
#   
#   res.ebs <- residuals(mod.ebs)
#   plot(fitted(mod.ebs), res.ebs, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals EBS")
#   
#   pp <- predict(mod.ebs, model.dat)
#   plot(model.dat$slp.sd, model.dat$ebs.sst.ar1, main = "Fitted EBS")
#   lines(model.dat$slp.sd, pp, col = "blue")
#   
#   mod.goa <- lm(goa.sst.ar1~slp.sd, data= model.dat)
#   goa.F <- glance(mod.goa)$statistic
#   goa.p <- glance(mod.goa)$p.value
#   
#   res.goa <- residuals(mod.goa)
#   plot(fitted(mod.goa), res.goa, xlab = "Fitted Values", ylab = "Residuals", main = "Residuals GOA")
#   
#   pp <- predict(mod.goa, model.dat)
#   plot(model.dat$slp.sd, model.dat$goa.sst.ar1, main = "Fitted GOA")
#   lines(model.dat$slp.sd, pp, col = "blue")
#   
#   obs.dat <- data.frame(region = c("EBS", "GOA"), 
#                         F_obs = c(ebs.F, goa.F), 
#                         p_obs = c(ebs.p, goa.p))
#   
#   sim.dat <- right_join(sim.dat, obs.dat) %>%
#               rename(F_sim = F_stat, p_sim = p_val)
#   
#   ## Calculate the number of times the simulations result in a greater F-stat and lower p than observed models
#   sim.dat %>%
#     group_by(region) %>%
#     reframe(larger_F = sum(F_sim > F_obs),
#             smaller_p = sum(p_sim < p_obs),
#             N = n(),
#             prop_largerF = larger_F/N,
#             prop_smallerp = smaller_p/N) %>%
#     dplyr::select(region, N, prop_largerF, prop_smallerp)
#     
#   # CESM2 ----
#   ## Read in CESM2 output files
#   sst <- rbind(readRDS(paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda")) %>%
#                  mutate(model = "FCM"),
#                readRDS(paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda")) %>%
#                  mutate(model = "MDM")) %>%
#     mutate(lat = as.numeric(lat), lon = as.numeric(lon))
#   
#   slp <- rbind(readRDS(paste0(dir, "Output/FCM_winterSLPa_ar1sd.rda")) %>%
#                  mutate(model = "FCM"),
#                readRDS(paste0(dir, "Output/MDM_winterSLPa_ar1sd.rda")) %>%
#                  mutate(model = "MDM")) %>%
#     mutate(lat = as.numeric(lat), lon = as.numeric(lon)) %>%
#     filter(member %in% sst$member)
#   
#   
#   ## Filter sst to EBS and GOA
#   # Specify EBS polygon
#   ebs.x <- c(183, 183, 203, 203, 191)
#   #ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
#   ebs.y <- c(53, 65, 65, 57.5, 53)
#   
#   # Filter sst initially to make things go quicker
#   sst.check <- sst %>%
#     filter(lon >= 170 & lon <= 210,
#            lat >= 50 & lat <= 70)
#   
#   # Create EBS polygon, filter sst not within that polygon
#   xp <- cbind(ebs.x, ebs.y)
#   loc= sst.check %>% dplyr::select(lon, lat)
#   check <- in.poly(loc, xp=xp)
#   
#   sst.ebs <- cbind(sst.check, check) %>%
#     filter(check == TRUE) %>%
#     dplyr::select(!check) 
#   
#   # Plot to check
#   sst.ebs %>%
#     filter(win.year == 1850, member == "1301.020", model == "FCM") -> tt
#   
#   ggplot()+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(50, 67), xlim = c(180, 205), expand = FALSE)+
#     geom_point(tt, mapping= aes(lon, lat))
#   
#   sst.ebs <- sst.ebs %>%
#     na.omit() %>%
#     group_by(win.year, member, model) %>%
#     reframe(ar1.winSSTa = mean(ar1.winSSTa))
#   
#   # Specify GOA polygon
#   goa.x <- c(201, 201, 205, 208, 225, 231, 201)
#   #goa.x <- ifelse(goa.x > 180, goa.x-360, goa.x)
#   goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)
#   
#   # Filter sst initially to make things go quicker
#   sst.check <- sst %>%
#     filter(lon >= 190 & lon <= 250,
#            lat >= 50 & lat <= 70)
#   
#   # Create EBS polygon, filter sst not within that polygon
#   xp <- cbind(goa.x, goa.y)
#   loc= sst.check %>% dplyr::select(lon, lat)
#   check <- in.poly(loc, xp=xp)
#   
#   sst.goa <- cbind(sst.check, check) %>%
#     filter(check == TRUE) %>%
#     dplyr::select(!check)
#   
#   # Plot to check
#   sst.goa %>%
#     filter(win.year == 1850, member == "1301.020", model == "FCM") -> tt
#   
#   ggplot()+
#     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
#     coord_cartesian(ylim = c(50, 67), xlim = c(190, 250), expand = FALSE)+
#     geom_point(tt, mapping= aes(lon, lat))
#   
#   sst.goa <- sst.goa %>%
#     na.omit() %>%
#     group_by(win.year, member, model) %>%
#     reframe(ar1.winSSTa = mean(ar1.winSSTa))
#   
#   
#   ## Get AR1 and SD of SLPa SD by ensemble member and model
#   SLP.pars <- slp %>%
#               na.omit() %>%
#               group_by(member, model) %>%
#               reframe(SLPaSD.ar1 = acf(sd.SLPa, plot = FALSE, na.action = na.pass)$acf[2],
#                       SLPaSD.sd = sd(sd.SLPa))
# 
#   
#   ## Create simulation dataframe
#   sim.dat <- data.frame()
#   
#   ## Create randomized SLPa SD timeseries, fit models for GOA and EBS SSTa
#   iter <- 1:10000
#   mem <- unique(SLP.pars$member)
#   sim.dat <- data.frame()
#  
#   for(ii in 1:length(iter)){
#     mod.dat <- data.frame()
#     
#     print(paste0("Iteration ", iter[ii], "/10000")) # progress tracking
#     
#     for(jj in 1:length(mem)){
#         # Filter params by member/model
#         SLP.pars %>%
#           filter(member == mem[jj]) -> sim.pars
#         
#         # Filter sst by member/model
#         sst.ebs %>%
#           filter(member == mem[jj]) -> sst.ebs2
#         
#         sst.goa %>%
#           filter(member == mem[jj]) -> sst.goa2
#         
#         # Simulate SLPa SD timeseries
#         SLPts <- arima.sim(model = list(order = c(1,0,0), 
#                                         ar = sim.pars$SLPaSD.ar1, 
#                                         sd= sim.pars$SLPaSD.sd),
#                            n = nrow(sst.ebs2)) # same # of rows for goa sst
#         
#         # Bind
#         out <- data.frame(member = sim.pars$member,
#                           model = sim.pars$model,
#                           win.year = sst.ebs2$win.year,
#                           ebs.sst.ar1 = sst.ebs2$ar1.winSSTa,
#                           goa.sst.ar1 = sst.goa2$ar1.winSSTa,
#                           slp.sd.sim = SLPts)
#         
#        mod.dat <- bind_rows(mod.dat, out)
#         
#      } # close ensemble member loop
#     
#     # Fit models
#     # EBS
#     mod.ebs <-lme(ebs.sst.ar1 ~ slp.sd.sim:model, data = mod.dat, random = ~ 1 | model/member)
#     
#     # GOA
#     mod.goa <-lme(goa.sst.ar1 ~ slp.sd.sim:model, data = mod.dat, random =  ~ 1 | model/member)
#     
#     # Create output df
#     out <- data.frame(iteration = iter[ii], 
#                       region = c("EBS", "EBS", "EBS", "GOA", "GOA", "GOA"), 
#                       param = rep(c("intercept", "slp.sd.sim:modelFCM", "slp.sd.sim:modelMDM"), 2),
#                       T_val = c(round(summary(mod.ebs)$tTable[,4],2), round(summary(mod.goa)$tTable[,4],2)),
#                       p_val = c(round(summary(mod.ebs)$tTable[,5],2), round(summary(mod.goa)$tTable[,5],2)))
#     
#     # Bind across iterations
#     sim.dat <- rbind(sim.dat, out)
#     
#   } # close iteration loop for simulations
#     
# 
#   write.csv(sim.dat, paste0(dir, "Output/CESM_lme_sim.csv"))
#   
#   sim.dat <- read.csv(paste0(dir, "Output/CESM_lme_sim.csv"))
#   
#   ## Fit models to observations and calculate T and p
#   slp2 <- slp %>%
#     na.omit() %>%
#     group_by(win.year, member, model) %>%
#     reframe(sd.winSLPa = mean(sd.SLPa))
#   
#   # EBS
#   right_join(slp2, sst.ebs) %>%
#     na.omit()-> mod.dat2
#   
#   mod.ebs <-lme(ar1.winSSTa ~ sd.winSLPa:model, data = mod.dat2, random = ~ 1 | model/member)
#   
#   # Effect plots
#   # compute partial residuals
#   mod.dat2$partial_resid <- residuals(mod.ebs) + predict(mod.ebs, level = 0)
#   
#   ggplot(mod.dat2, aes(x = sd.winSLPa, y = partial_resid)) +
#     geom_point(alpha = 0.6) + 
#     geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +  # Fixed-effect slope
#     facet_wrap(~ member) +  # Facet by member
#     ggtitle("EBS")+
#     labs(x = "sd.winSLPa", y = "Partial Residuals") +
#     theme_bw() -> ebs.resid
#   
#   ggsave(plot = ebs.resid, "./Figures/ebsCESM.partialresid.png", width = 8, height = 8, units= 'in')
#   
#   
#   
#   rr <- Effect(c("sd.winSLPa", "model"), mod.ebs, partial.residuals = TRUE)
#   plot(rr,  main = "EBS partial residuals", show_residuals = TRUE, show_residuals_line = FALSE)
#   
#   rr <- Effect(c("sd.winSLPa", "model"), mod.ebs, confidence.level = 0.95)
#   plot(rr, main = "EBS fitted values with 95% CI")
#   
#   tt <- predictorEffect("sd.winSLPa", mod.ebs)
#   
#   # GOA
#   right_join(slp2, sst.goa) %>%
#     na.omit()-> mod.dat2
#   
#   mod.goa1 <-lme(ar1.winSSTa ~ sd.winSLPa:model, data = mod.dat2, random = ~ 1 | model/member)
#   mod.goa2 <-lme(ar1.winSSTa ~ sd.winSLPa:model, weights = varIdent(), data = mod.dat2, random = ~ 1 | model/member)
#   mod.goa <-lme(ar1.winSSTa ~ sd.winSLPa*model, data = mod.dat2, random = ~ 1 | member)
#   
#   
#   
#   # Effect plots
#   # compute partial residuals
#   mod.dat2$partial_resid <- residuals(mod.goa) + predict(mod.goa, level = 0)
#   
#   ggplot(mod.dat2, aes(x = sd.winSLPa, y = partial_resid)) +
#     geom_point(alpha = 0.6) + 
#     geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +  # Fixed-effect slope
#     facet_wrap(~ member) +  # Facet by member
#     geom_smooth(method = "loess", se = FALSE) +  # Smoothed residuals
#     ggtitle("GOA")+
#     labs(x = "sd.winSLPa", y = "Partial Residuals") + #(Fixed Effect + Residuals)
#     theme_bw() -> goa.resid
#   
#   ggsave(plot = goa.resid, "./Figures/goaCESM.partialresid.png", width = 8, height = 8, units= 'in')
#   
#   rr <- Effect(c("sd.winSLPa", "model"), mod.goa, partial.residuals = TRUE)
#   plot(rr,  main = "GOA partial residuals", show_residuals = TRUE, show_residuals_line = FALSE)
#   
#   rr <- Effect(c("sd.winSLPa", "model"), mod.goa, confidence.level = 0.95)
#   plot(rr, main = "GOA fitted values with 95% CI")
#   
#   tt <- predictorEffect("sd.winSLPa", mod.goa)
#   
#   # Pull out observed model params
#   obs.dat <- data.frame(region = c("EBS", "EBS", "EBS", "GOA", "GOA", "GOA"), 
#                     param = rep(c("intercept", "slp.sd.sim:modelFCM", "slp.sd.sim:modelMDM"), 2),
#                     T_obs = c(round(summary(mod.ebs)$tTable[,4],2), round(summary(mod.goa)$tTable[,4],2)),
#                     p_obs = c(round(summary(mod.ebs)$tTable[,5],2), round(summary(mod.goa)$tTable[,5],2)))
#   
#   sim.dat2 <- right_join(sim.dat, obs.dat) %>%
#     rename(T_sim = T_val, p_sim = p_val)
#   
#   
#   ## Calculate the number of times the simulations result in a greater F-stat and lower p than observed models
#   sim.dat2 %>%
#     group_by(region, param) %>%
#     reframe(larger_T = sum(T_sim > T_obs),
#             smaller_p = sum(p_sim < p_obs),
#             N = n(),
#             prop_largerT = larger_T/N,
#             prop_smallerp = smaller_p/N) %>%
#     dplyr::select(region, param, N, prop_largerT, prop_smallerp) %>%
#     filter(param!="intercept") %>%
#     rename(term = param)
#   
#   ## Compare different model parameterizations
#   # EBS
#   right_join(slp2, sst.ebs) %>%
#     na.omit()-> mod.dat2
#   
#   mod.ebs1 <-lme(ar1.winSSTa ~ sd.winSLPa:model, data = mod.dat2, random = ~ 1 | model/member, method = "REML") #
#   mod.ebs2 <-lme(ar1.winSSTa ~ sd.winSLPa*model, data = mod.dat2, random = ~ 1 | member, method = "REML") # intercept by model
#   
#   AICc(mod.ebs1, mod.ebs2)
#   
#   
#   right_join(slp2, sst.goa) %>%
#     na.omit()-> mod.dat2
#   
#   mod.goa1 <-lme(ar1.winSSTa ~ sd.winSLPa:model, data = mod.dat2, random = ~ 1 | model/member, method = "REML") #
#   mod.goa2 <-lme(ar1.winSSTa ~ sd.winSLPa*model, data = mod.dat2, random = ~ 1 | member, method = "REML") # intercept by model
#   
#   AICc(mod.goa1, mod.goa2)
#   
#   mod.ebs2 <- lme(
#               ar1.winSSTa ~ sd.winSLPa * model,  # Interaction term (model-specific slopes)
#               random = ~1 | model/member,
#               data = mod.dat2,
#               method = "ML"
#             )
#   
#   summary(mod.ebs2)
#   
# 7) CESM EOFs by ensemble member ----
  # Identify folders to download files from
  fcm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/FCM/") #FCM SST
  mdm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/MDM/") #MDM SST
  
  fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
  mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP
  
  # Extract time info for processing below (same across files)
  files <- list.files(fcm.sst.dir, full.names = TRUE)
  
  time_units <- ncmeta::nc_atts(files[1], "time") %>% 
    filter(name == "units") %>%
    pull(value)
  
  unit_parts <- str_split(time_units, " since ")[[1]]
  time_unit <- unit_parts[1]
  origin_date <- ymd_hms(unit_parts[2])
  
  # FCM SST ----  
  files <- list.files(fcm.sst.dir, full.names = TRUE)
  
    # FIRST CHUNK
    files <- files[c(1:25)]
    
    fcm.sst.EOF1.1 <- data.table()
    fcm.sst.PC1.1 <- data.frame()
    
    #files <- files[1:2]
    
    for(ii in 1:length(files)){
      gc()
      print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
      
      # load and process file
      tidync(files[ii]) %>%
        hyper_filter(lon = lon >= 125 & lon <= 255,
                     lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
        #time = time > 711475) %>% # greater than 1947
        activate("SST") %>%
        hyper_tibble() %>%
        mutate(time = origin_date + lubridate::days(time),
               year = lubridate::year(time),
               month = lubridate::month(time),
               member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
        group_by(lat, lon, month, member) %>%
        mutate(mean.month.SST = mean(SST),
               sd.month.SST = sd(SST)) %>% # compute monthly mean by grid cell and member
        ungroup() %>%
        mutate(SSTa.deseasoned = SST-mean.month.SST, # remove monthly climatology
               time.index = row_number()) -> out # compute anomalies/z-scores
      
      
      # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
      setDT(out) # convert to data.table
      
      out[, detrended.SSTa := { # output data table column
        fit <- lm(SSTa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
        residuals(fit)           # Extract residuals as detrended values
      }, by = .(lat, lon, member)]
      
      
      # Scale each grid point by its SD
      out[, detrended.SSTa.scaled := { # output data table column
        scale(detrended.SSTa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
      }, by = .(lat, lon, member)]
      
      
     # Format data for EOFs
      out %>%
        dplyr::select(time, lon, lat, member, detrended.SSTa.scaled) %>%
        distinct(.) %>%
        as.data.frame(.) %>%
        unite("crds", lat, lon, sep = "-") %>%
        pivot_wider(., names_from= crds, values_from = detrended.SSTa.scaled) -> EOF.dat
      
    # Calculate EOF weights
      lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
      
      length(lat)
      ncol(EOF.dat[,-c(1:2)])
      
      weights <- sqrt(cos(lat*pi/180))
      
    # Run EOF
      print("Running EOF")
      eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
      
  
    # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
      # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
      U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
      
      PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
             as.data.frame() %>%
             mutate(time = EOF.dat$time,
                    V1 = scale(as.numeric(V1)),
                    member = unique(EOF.dat$member)) %>%
        rename(PC1 = V1) %>%
        mutate(PC1= as.numeric(PC1))
      
      #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
      
      V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
      
      #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
      
      # output dataframe
      EOF1.dat <- expand.grid(
        location = colnames(EOF.dat[, -c(1,2)]),
        stringsAsFactors = FALSE) %>%
        mutate(member= unique(EOF.dat$member),
               EOF1 = V.1[match(.$location, colnames(EOF.dat))])
        
      setDT(EOF1.dat)

      fcm.sst.EOF1.1 <- bind_rows(fcm.sst.EOF1.1, EOF1.dat)
      fcm.sst.PC1.1 <- bind_rows(fcm.sst.PC1.1, PC1.dat)
      
    }
    
    # SECOND CHUNK
    files <- list.files(fcm.sst.dir, full.names = TRUE)
    files <- files[c(26:48)] # last two files result in empty slice
    
    fcm.sst.EOF1.2 <- data.table()
    fcm.sst.PC1.2 <- data.table()
    
    for(ii in 1:length(files)){
      gc()
      print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
      
      # load and process file
      tidync(files[ii]) %>%
        hyper_filter(lon = lon >= 125 & lon <= 255,
                     lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
        #time = time > 711475) %>% # greater than 1947
        activate("SST") %>%
        hyper_tibble() %>%
        mutate(time = origin_date + lubridate::days(time),
               year = lubridate::year(time),
               month = lubridate::month(time),
               member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
        group_by(lat, lon, month, member) %>%
        mutate(mean.month.SST = mean(SST),
               sd.month.SST = sd(SST)) %>% # compute monthly mean by grid cell and member
        ungroup() %>%
        mutate(SSTa.deseasoned = SST-mean.month.SST, # remove monthly climatology
               time.index = row_number()) -> out # compute anomalies/z-scores
      
      
      # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
      setDT(out) # convert to data.table
      
      out[, detrended.SSTa := { # output data table column
        fit <- lm(SSTa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
        residuals(fit)           # Extract residuals as detrended values
      }, by = .(lat, lon, member)]
      
      
      # Scale each grid point by its SD
      out[, detrended.SSTa.scaled := { # output data table column
        scale(detrended.SSTa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
      }, by = .(lat, lon, member)]
      
      
      # Format data for EOFs
      out %>%
        dplyr::select(time, lon, lat, member, detrended.SSTa.scaled) %>%
        distinct(.) %>%
        as.data.frame(.) %>%
        unite("crds", lat, lon, sep = "-") %>%
        pivot_wider(., names_from= crds, values_from = detrended.SSTa.scaled) -> EOF.dat
      
      # Calculate EOF weights
      lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
      
      length(lat)
      ncol(EOF.dat[,-c(1:2)])
      
      weights <- sqrt(cos(lat*pi/180))
      
      # Run EOF
      print("Running EOF")
      eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
      
      
      # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
      # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
      U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
      
      PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
        as.data.frame() %>%
        mutate(time = EOF.dat$time,
               V1 = scale(as.numeric(V1)),
               member = unique(EOF.dat$member)) %>%
        rename(PC1 = V1) %>%
        mutate(PC1= as.numeric(PC1))
      
      #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
      
      V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
      
      #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
      
      # output dataframe
      EOF1.dat <- expand.grid(
        location = colnames(EOF.dat[, -c(1,2)]),
        stringsAsFactors = FALSE) %>%
        mutate(member= unique(EOF.dat$member),
               EOF1 = V.1[match(.$location, colnames(EOF.dat))])
      
      setDT(EOF1.dat)
      
      fcm.sst.EOF1.2 <- bind_rows(fcm.sst.EOF1.2, EOF1.dat)
      fcm.sst.PC1.2<- bind_rows(fcm.sst.PC1.2, PC1.dat)
      
    }
    
    # BIND CHUNKS
    rbind(fcm.sst.EOF1.1, fcm.sst.EOF1.2) -> fcm.sst.EOF1
    rbind(fcm.sst.PC1.1, fcm.sst.PC1.2) -> fcm.sst.PC1
    
    setDT(fcm.sst.EOF1)
    setDT(fcm.sst.PC1)
    
    # Save file
    saveRDS(fcm.sst.EOF1, paste0(dir, "Output/FCM_SSTa_EOF1.rda"))
    saveRDS(fcm.sst.PC1, paste0(dir, "Output/FCM_SSTa_PC1.rda"))
    
    
  # MDM SST ----
  files <- list.files(mdm.sst.dir, full.names = TRUE)
  mdm.sst.EOF1 <- tibble()
  mdm.sst.PC1 <- tibble()
  
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
      #time = time > 711475) %>% # greater than 1947
      activate("SST") %>%
      hyper_tibble() %>%
      mutate(time = origin_date + lubridate::days(time),
             year = lubridate::year(time),
             month = lubridate::month(time),
             member = substr(files[ii], 68, 78),
             member= gsub(".postP.", ".", member)) %>% # extracting ensemble member #
      mutate(mean.month.SST = mean(SST),
             sd.month.SST = sd(SST)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SSTa.deseasoned = SST-mean.month.SST, # remove monthly climatology
             time.index = row_number()) -> out # compute anomalies/z-scores
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended.SSTa := { # output data table column
      fit <- lm(SSTa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]
    
    
    # Scale each grid point by its SD
    out[, detrended.SSTa.scaled := { # output data table column
      scale(detrended.SSTa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
    }, by = .(lat, lon, member)]
    
    
    # Format data for EOFs
    out %>%
      dplyr::select(time, lon, lat, member, detrended.SSTa.scaled) %>%
      distinct(.) %>%
      as.data.frame(.) %>%
      unite("crds", lat, lon, sep = "-") %>%
      pivot_wider(., names_from= crds, values_from = detrended.SSTa.scaled) -> EOF.dat
    
    # Calculate EOF weights
    lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
    
    length(lat)
    ncol(EOF.dat[,-c(1:2)])
    
    weights <- sqrt(cos(lat*pi/180))
    
    # Run EOF
    print("Running EOF")
    eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
    
    
    # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
    # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
    U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
    
    PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
      as.data.frame() %>%
      mutate(time = EOF.dat$time,
             V1 = scale(as.numeric(V1)),
             member = unique(EOF.dat$member)) %>%
      rename(PC1 = V1) %>%
      mutate(PC1= as.numeric(PC1))
    
    #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
    
    V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
    
    #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
    
    # output dataframe
    EOF1.dat <- expand.grid(
      location = colnames(EOF.dat[, -c(1,2)]),
      stringsAsFactors = FALSE) %>%
      mutate(member= unique(EOF.dat$member),
             EOF1 = V.1[match(.$location, colnames(EOF.dat))])
    
    setDT(EOF1.dat)
    
    mdm.sst.EOF1 <- bind_rows(mdm.sst.EOF1, EOF1.dat)
    mdm.sst.PC1 <- bind_rows(mdm.sst.PC1, PC1.dat)
  }
  
  setDT(mdm.sst.EOF1)
  setDT(mdm.sst.PC1)
  
  
  # Save file
  saveRDS(mdm.sst.EOF1, paste0(dir, "Output/MDM_SSTa_EOF1.rda"))
  saveRDS(mdm.sst.PC1, paste0(dir, "Output/MDM_SSTa_PC1.rda"))
  
  
  # FCM SLP ----
    # FIRST CHUNK
    files <- list.files(fcm.slp.dir, full.names = TRUE)
   
    files <- files[c(1:25)]
  
    fcm.slp.EOF1.1 <- tibble()
    fcm.slp.PC1.1 <- tibble()
    
    
    for(ii in 1:length(files)){
      
      print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
      
      # load and process file
      tidync(files[ii]) %>%
        hyper_filter(lon = lon >= 125 & lon <= 255,
                     lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
        #time = time > 711475) %>% # greater than 1947
        activate("PSL") %>%
        hyper_tibble() %>%
        mutate(year = lubridate::year(time), # already in correct format
               month = lubridate::month(time),
               day = lubridate::day(time),
               member = substr(files[ii], 78, 85)) %>% # extracting ensemble member #
        mutate(mean.month.SLP = mean(PSL),
               sd.month.SLP = sd(PSL)) %>% # compute monthly mean by grid cell and member
        ungroup() %>%
        mutate(SLPa.deseasoned = PSL-mean.month.SLP, # remove monthly climatology
               time.index = row_number()) -> out # compute anomalies/z-scores
      
      
      # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
      setDT(out) # convert to data.table
      
      out[, detrended.SLPa := { # output data table column
        fit <- lm(SLPa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
        residuals(fit)           # Extract residuals as detrended values
      }, by = .(lat, lon, member)]
      
      
      # Scale each grid point by its SD
      out[, detrended.SLPa.scaled := { # output data table column
        scale(detrended.SLPa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
      }, by = .(lat, lon, member)]
      
      
      # Format data for EOFs
      out %>%
        dplyr::select(time, lon, lat, member, detrended.SLPa.scaled) %>%
        distinct(.) %>%
        as.data.frame(.) %>%
        unite("crds", lat, lon, sep = "-") %>%
        pivot_wider(., names_from= crds, values_from = detrended.SLPa.scaled) -> EOF.dat
      
      # Calculate EOF weights
      lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
      
      length(lat)
      ncol(EOF.dat[,-c(1:2)])
      
      weights <- sqrt(cos(lat*pi/180))
      
      # Run EOF
      print("Running EOF")
      eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
      
      
      # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
      # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
      U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
      
      PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
        as.data.frame() %>%
        mutate(time = EOF.dat$time,
               V1 = scale(as.numeric(V1)),
               member = unique(EOF.dat$member)) %>%
        rename(PC1 = V1) %>%
        mutate(PC1= as.numeric(PC1))
      
      #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
      
      V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
      
      #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
      
      # output dataframe
      EOF1.dat <- expand.grid(
        location = colnames(EOF.dat[, -c(1,2)]),
        stringsAsFactors = FALSE) %>%
        mutate(member= unique(EOF.dat$member),
               EOF1 = V.1[match(.$location, colnames(EOF.dat))])
      
      setDT(EOF1.dat)
      
      fcm.slp.EOF1.1 <- bind_rows(fcm.slp.EOF1.1, EOF1.dat)
      fcm.slp.PC1.1<- bind_rows(fcm.slp.PC1.1, PC1.dat)
      
    }
    
    # SECOND CHUNK
    files <- list.files(fcm.slp.dir, full.names = TRUE)
    
    files <- files[c(26:48)] # last two files result in empty slice
    
    fcm.slp.EOF1.2 <- tibble()
    fcm.slp.PC1.2 <- tibble()
    
    
    for(ii in 1:length(files)){
      
      print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
      
      # load and process file
      tidync(files[ii]) %>%
        hyper_filter(lon = lon >= 125 & lon <= 255,
                     lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
        #time = time > 711475) %>% # greater than 1947
        activate("PSL") %>%
        hyper_tibble() %>%
        mutate(year = lubridate::year(time), # already in correct format
               month = lubridate::month(time),
               day = lubridate::day(time),
               member = substr(files[ii], 78, 85)) %>% # extracting ensemble member #
        mutate(mean.month.SLP = mean(PSL),
               sd.month.SLP = sd(PSL)) %>% # compute monthly mean by grid cell and member
        ungroup() %>%
        mutate(SLPa.deseasoned = PSL-mean.month.SLP, # remove monthly climatology
               time.index = row_number()) -> out # compute anomalies/z-scores
      
      
      # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
      setDT(out) # convert to data.table
      
      out[, detrended.SLPa := { # output data table column
        fit <- lm(SLPa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
        residuals(fit)           # Extract residuals as detrended values
      }, by = .(lat, lon, member)]
      
      
      # Scale each grid point by its SD
      out[, detrended.SLPa.scaled := { # output data table column
        scale(detrended.SLPa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
      }, by = .(lat, lon, member)]
      
      
      # Format data for EOFs
      out %>%
        dplyr::select(time, lon, lat, member, detrended.SLPa.scaled) %>%
        distinct(.) %>%
        as.data.frame(.) %>%
        unite("crds", lat, lon, sep = "-") %>%
        pivot_wider(., names_from= crds, values_from = detrended.SLPa.scaled) -> EOF.dat
      
      # Calculate EOF weights
      lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
      
      length(lat)
      ncol(EOF.dat[,-c(1:2)])
      
      weights <- sqrt(cos(lat*pi/180))
      
      # Run EOF
      print("Running EOF")
      eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
      
      
      # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
      # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
      U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
      
      PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
        as.data.frame() %>%
        mutate(time = EOF.dat$time,
               V1 = scale(as.numeric(V1)),
               member = unique(EOF.dat$member)) %>%
        rename(PC1 = V1) %>%
        mutate(PC1= as.numeric(PC1))
      
      #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
      
      V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
      
      #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
      
      # output dataframe
      EOF1.dat <- expand.grid(
        location = colnames(EOF.dat[, -c(1,2)]),
        stringsAsFactors = FALSE) %>%
        mutate(member= unique(EOF.dat$member),
               EOF1 = V.1[match(.$location, colnames(EOF.dat))])
      
      setDT(EOF1.dat)
      
      fcm.slp.EOF1.2 <- bind_rows(fcm.slp.EOF1.2, EOF1.dat)
      fcm.slp.PC1.2<- bind_rows(fcm.slp.PC1.2, PC1.dat)
      
    }
    
    # BIND CHUNKS
    rbind(fcm.slp.EOF1.1, fcm.slp.EOF1.2) -> fcm.slp.EOF1
    rbind(fcm.slp.PC1.1, fcm.slp.PC1.2) -> fcm.slp.PC1
    
    setDT(fcm.slp.EOF1)
    setDT(fcm.slp.PC1)
    
    # Save file
    saveRDS(fcm.slp.EOF1, paste0(dir, "Output/FCM_slpa_EOF1.rda"))
    saveRDS(fcm.slp.PC1, paste0(dir, "Output/FCM_slpa_PC1.rda"))
    
  
  # MDM SLP ----
  files <- list.files(mdm.slp.dir, full.names = TRUE)
  mdm.slp.EOF <- tibble()
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 20 & lat <= 68) %>% # extra tropical north pacific region
      #time = time > 711475) %>% # greater than 1947
      activate("PSL") %>%
      hyper_tibble() %>%
      mutate(year = lubridate::year(time), # already in correct format
             month = lubridate::month(time),
             day = lubridate::day(time),
             member = substr(files[ii], 68, 72)) %>% # extracting ensemble member #
      mutate(mean.month.SLP = mean(PSL),
             sd.month.SLP = sd(PSL)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SLPa.deseasoned = PSL-mean.month.SLP, # remove monthly climatology
             time.index = row_number()) -> out # compute anomalies/z-scores
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended.SLPa := { # output data table column
      fit <- lm(SLPa.deseasoned ~ time.index)  # Fit linear model for each lat/lon group across months/years
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon, member)]
    
    
    # Scale each grid point by its SD
    out[, detrended.SLPa.scaled := { # output data table column
      scale(detrended.SLPa, center = FALSE, scale = TRUE) # center = FALSE bc data is deseasonlized and detrended (mean close to zero); scale = TRUE divides by standard deviation
    }, by = .(lat, lon, member)]
    
    
    # Format data for EOFs
    out %>%
      dplyr::select(time, lon, lat, member, detrended.SLPa.scaled) %>%
      distinct(.) %>%
      as.data.frame(.) %>%
      unite("crds", lat, lon, sep = "-") %>%
      pivot_wider(., names_from= crds, values_from = detrended.SLPa.scaled) -> EOF.dat
    
    # Calculate EOF weights
    lat <- as.numeric(sapply(strsplit(colnames(EOF.dat[,-c(1:2)]), "-"), `[`, 1))
    
    length(lat)
    ncol(EOF.dat[,-c(1:2)])
    
    weights <- sqrt(cos(lat*pi/180))
    
    # Run EOF
    print("Running EOF")
    eof <- FactoMineR::svd.triplet(cov(EOF.dat[,-c(1:2)]), col.w=na.omit(weights)) #weighting the columns
    
    
    # Pull out loadings and pc scores, scale, bind to lat/lon. Because our rows represent time and columns represent location in 
    # EOF.dat, U represents temporal patterns (PC scores) and V represents spatial patterns (loadings)
    U.1 <- scale(eof$U[,1] * eof$vs[1]) # temporal pattern (PC1 score) 
    
    PC1.dat <- as.matrix(EOF.dat[, -c(1:2)]) %*% U.1 %>%
      as.data.frame() %>%
      mutate(time = EOF.dat$time,
             V1 = scale(as.numeric(V1)),
             member = unique(EOF.dat$member)) %>%
      rename(PC1 = V1) %>%
      mutate(PC1= as.numeric(PC1))
    
    #U.2 <- scale(eof$U[,2] * eof$vs[2]) # temporal pattern (PC2 score)
    
    V.1 <- scale(eof$V[,1]) # spatial pattern (EOF1 loading)
    
    #V.2 <- scale(eof$V[,2]) # spatial pattern (EOF2 loading)
    
    # output dataframe
    EOF1.dat <- expand.grid(
      location = colnames(EOF.dat[, -c(1,2)]),
      stringsAsFactors = FALSE) %>%
      mutate(member= unique(EOF.dat$member),
             EOF1 = V.1[match(.$location, colnames(EOF.dat))])
    
    setDT(EOF1.dat)
    
    mdm.slp.EOF1 <- bind_rows(mdm.slp.EOF1, EOF1.dat)
    mdm.slp.PC1<- bind_rows(mdm.slp.PC1, PC1.dat)
    
  }
  
  setDT(mdm.slp.EOF1)
  setDT(mdm.slp.PC1)
  
  
  # Save file
  saveRDS(mdm.slp.EOF1, paste0(dir, "Output/MDM_SLPa_EOF1.rda"))
  saveRDS(mdm.slp.PC1, paste0(dir, "Output/MDM_SLPa_PC1.rda"))
  
  
  
  # PLOT EOF1 loadings and PC1 scores----
  # fcm sst
    fcm.sst.eof <- readRDS(paste0(dir, "Output/FCM_SSTa_EOF1.rda")) %>%
                    mutate(lat = as.numeric(sapply(strsplit(location, "-"), `[`, 1)),
                           lon = as.numeric(sapply(strsplit(location, "-"), `[`, 2)))
    fcm.sst.pc <- readRDS(paste0(dir, "Output/FCM_SSTa_PC1.rda"))
    
    fcm.sst.pc %>%
      mutate(year = year(time)) -> pp
    
    # EOF1
    ggplot()+
     geom_tile(fcm.sst.eof, mapping= aes(lon, lat, fill = EOF1))+
     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
     coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
     xlab("Latitude")+
     facet_wrap(~member)+
     ggtitle("FCM SSTa EOF1")+
     ylab("Longitude")+
     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                          midpoint=  median(na.omit(fcm.sst.eof$EOF1)),
                          name = "EOF 1",
                          limits = c(-(max(na.omit(abs(fcm.sst.eof$EOF1)))), max(na.omit(abs(fcm.sst.eof$EOF1)))))
   
   ggsave("./Figures/FCM.SSTa.EOF1.png", width = 11, height=  8.5, units = "in")
   
   # PC1
   ggplot()+
     geom_line(fcm.sst.pc, mapping= aes(time, PC1), linewidth = 0.15)+
     xlab("Time")+
     facet_wrap(~member)+
     ylab("PC1")+
     theme_bw()+
     ggtitle("FCM SSTa PC1")+
     theme(axis.text  = element_text(size = 10),
           axis.title = element_text(size = 11))
   
   ggsave("./Figures/FCM.SSTa.PC1.png", width = 11, height=  8.5, units = "in")
   

 # mdm sst
   mdm.sst.eof <- readRDS(paste0(dir, "Output/MDM_SSTa_EOF1.rda")) %>%
     mutate(lat = as.numeric(sapply(strsplit(location, "-"), `[`, 1)),
            lon = as.numeric(sapply(strsplit(location, "-"), `[`, 2)))
   mdm.sst.pc <- readRDS(paste0(dir, "Output/MDM_SSTa_PC1.rda"))
   
   # EOF1
   ggplot()+
     geom_tile(mdm.sst.eof, mapping= aes(lon, lat, fill = EOF1))+
     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
     coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
     xlab("Latitude")+
     facet_wrap(~member)+
     ggtitle("MDM SSTa EOF1")+
     ylab("Longitude")+
     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                          midpoint=  median(na.omit(mdm.sst.eof$EOF1)),
                          name = "EOF 1",
                          limits = c(-(max(na.omit(abs(mdm.sst.eof$EOF1)))), max(na.omit(abs(mdm.sst.eof$EOF1)))))
   
   ggsave("./Figures/MDM.SSTa.EOF1.png", width = 11, height=  8.5, units = "in")
   
   # PC1
   ggplot()+
     geom_line(mdm.sst.pc, mapping= aes(time, PC1), linewidth = 0.15)+
     xlab("Time")+
     facet_wrap(~member)+
     ylab("PC1")+
     theme_bw()+
     ggtitle("MDM SSTa PC1")+
     theme(axis.text  = element_text(size = 10),
           axis.title = element_text(size = 11))
   
   ggsave("./Figures/MDM.SSTa.PC1.png", width = 11, height=  8.5, units = "in")
   
  
  # fcm slp
   fcm.slp.eof <- readRDS(paste0(dir, "Output/FCM_SLPa_EOF1.rda")) %>%
     mutate(lat = as.numeric(sapply(strsplit(location, "-"), `[`, 1)),
            lon = as.numeric(sapply(strsplit(location, "-"), `[`, 2)))
   fcm.slp.pc <- readRDS(paste0(dir, "Output/FCM_SLPa_PC1.rda"))
   
   fcm.slp.pc %>%
     mutate(year = year(time)) -> pp
   
   # EOF1
   ggplot()+
     geom_tile(fcm.slp.eof, mapping= aes(lon, lat, fill = EOF1))+
     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
     coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
     xlab("Latitude")+
     facet_wrap(~member)+
     ggtitle("FCM SLPa EOF1")+
     ylab("Longitude")+
     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                          midpoint=  median(na.omit(fcm.slp.eof$EOF1)),
                          name = "EOF 1",
                          limits = c(-(max(na.omit(abs(fcm.slp.eof$EOF1)))), max(na.omit(abs(fcm.slp.eof$EOF1)))))
   
   ggsave("./Figures/FCM.SLPa.EOF1.png", width = 11, height=  8.5, units = "in")
   
   # PC1
   ggplot()+
     geom_line(fcm.slp.pc, mapping= aes(time, PC1), linewidth = 0.15)+
     xlab("Time")+
     facet_wrap(~member)+
     ylab("PC1")+
     theme_bw()+
     ggtitle("FCM SLPa PC1")+
     theme(axis.text  = element_text(size = 10),
           axis.title = element_text(size = 11))
   
   ggsave("./Figures/FCM.SLPa.PC1.png", width = 11, height=  8.5, units = "in")
   
   # mdm slp
   mdm.slp.eof <- readRDS(paste0(dir, "Output/MDM_SLPa_EOF1.rda")) %>%
     mutate(lat = as.numeric(sapply(strsplit(location, "-"), `[`, 1)),
            lon = as.numeric(sapply(strsplit(location, "-"), `[`, 2)))
   mdm.slp.pc <- readRDS(paste0(dir, "Output/MDM_SLPa_PC1.rda"))
   
   mdm.slp.pc %>%
     mutate(year = year(time)) -> pp
   
   # EOF1
   ggplot()+
     geom_tile(mdm.slp.eof, mapping= aes(lon, lat, fill = EOF1))+
     geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
     coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
     xlab("Latitude")+
     facet_wrap(~member)+
     ggtitle("MDM SLPa EOF1")+
     ylab("Longitude")+
     scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white",
                          midpoint=  median(na.omit(mdm.slp.eof$EOF1)),
                          name = "EOF 1",
                          limits = c(-(max(na.omit(abs(mdm.slp.eof$EOF1)))), max(na.omit(abs(mdm.slp.eof$EOF1)))))
   
   ggsave("./Figures/MDM.SLPa.EOF1.png", width = 11, height=  8.5, units = "in")
   
   # PC1
   ggplot()+
     geom_line(mdm.slp.pc, mapping= aes(time, PC1), linewidth = 0.15)+
     xlab("Time")+
     facet_wrap(~member)+
     ylab("PC1")+
     theme_bw()+
     ggtitle("MDM SLPa PC1")+
     theme(axis.text  = element_text(size = 10),
           axis.title = element_text(size = 11))
   
   ggsave("./Figures/MDM.SLPa.PC1.png", width = 11, height=  8.5, units = "in")
   
   