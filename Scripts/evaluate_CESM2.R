# PURPOSE: to download CESM2 fcm and mdm SST/SLP model outputs and calculate ar1/SD

# Note: steps 1) and 2) take some time to process. If you just need to modify the final outputs/figures, proceed
# to step 3)

# Author: Emily Ryznar

# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# Get map layers
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))


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
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking
    
    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 0 & lat <= 75) %>% # north pacific region
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
      mutate(SSTa = SST - mean.month.SST) %>% # compute anomalies
      group_by(lon, lat, year, member) %>% 
      reframe(mean.SSTa = mean(SSTa))-> out # calculate mean annual SSTa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table

    out[, detrended_SSTa := { # output data table column
      fit <- lm(mean.SSTa ~ year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon)]
    
    # Calculate grid cell AR1
    out %>%
      group_by(lon, lat, member) %>%
      reframe(ar1.SSTa = sapply(rollapply(detrended_SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2), # ar1 over 15 year rolling windows
              sd.SSTa = rollapply(detrended_SSTa, width = 15, FUN = sd), # sd over 15-year rolling windows
              ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
              ar1.mean = mean(ar1.SSTa), # mean AR1 of rolling window timeseries
              sd.sd = sd(sd.SSTa), # sd of SD of rolling window timeseries
              sd.mean = mean(sd.SSTa)) %>% # mean SD of rolling window timeseries
      mutate(lon = as.numeric(lon),
             lat = as.numeric(lat)) %>%
      dplyr::select(!c(ar1.SSTa, sd.SSTa)) %>%
      distinct() -> ar1.sd.dat
    
    # stack processed files
    fcm.sst <- bind_rows(fcm.sst, ar1.sd.dat)
    
  }

  # Save file
  saveRDS(fcm.sst, paste0(dir, "Output/processed.fcm.sst.rda"))
  
  # MDM SST ----
  files <- list.files(mdm.sst.dir, full.names = TRUE)
  mdm.sst <- tibble()
  
  for(ii in 1:length(files)){
    
    print(paste0("Processing file ", (1:length(files))[ii], "/", length(files))) # for progress tracking

    # load and process file
    tidync(files[ii]) %>%
      hyper_filter(lon = lon >= 125 & lon <= 255,
                   lat = lat >= 0 & lat <= 75) %>% # north pacific region
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
      mutate(SSTa = SST - mean.month.SST) %>% # compute anomalies
      group_by(lon, lat, year, member) %>% 
      reframe(mean.SSTa = mean(SSTa)) -> out # calculate mean annual SSTa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SSTa := { # output column name
      fit <- lm(mean.SSTa ~ year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon)]
    
    # Calculate grid cell AR1
    out %>%
      group_by(lon, lat, member) %>%
      reframe(ar1.SSTa = sapply(rollapply(detrended_SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2), # ar1 over 15 year rolling windows
              sd.SSTa = rollapply(detrended_SSTa, width = 15, FUN = sd), # sd over 15-year rolling windows
              ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
              ar1.mean = mean(ar1.SSTa), # mean AR1 of rolling window timeseries
              sd.sd = sd(sd.SSTa), # sd of SD of rolling window timeseries
              sd.mean = mean(sd.SSTa)) %>% # mean SD of rolling window timeseries
      mutate(lon = as.numeric(lon),
             lat = as.numeric(lat)) %>%
      dplyr::select(!c(ar1.SSTa, sd.SSTa)) %>%
      distinct() -> ar1.sd.dat
    
    # stack processed files
    mdm.sst <- bind_rows(mdm.sst, ar1.sd.dat)
    
    
  }
  
  # Save file
  saveRDS(mdm.sst, paste0(dir, "Output/processed.mdm.sst.rda"))
  
  
  
# 3) CALCULATE/PLOT CELL-WISE AR1 SD, SST SD, and MEAN AR1 ACROSS ENSEMBLE ------------------------
# Calculate mean SD AR1 and AR1 by grid cell across ensemble members
fcm.sst <- readRDS(paste0(dir, "Output/processed.fcm.sst.rda"))

fcm.sst %>%
  filter(lat >=20 & lat<=68) %>% # isolate extra tropical N. pacific
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Fully coupled") -> fcm.sst.ar1sd

mdm.sst <- readRDS(paste0(dir, "Output/processed.mdm.sst.rda"))

mdm.sst %>%
  filter(lat >=20 & lat<=68) %>% # isolate extra tropical N. pacific
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Mechanically decoupled") -> mdm.sst.ar1sd

# Join
plot.dat <- rbind(fcm.sst.ar1sd, mdm.sst.ar1sd)

# Calculate differences
plot.dat %>%
  # group_by(type) %>%
  # mutate(ar1.sd = scale(ar1.sd)[,1],
  #        ar1.mean = scale(ar1.mean)[,1],
  #        sd.mean = scale(sd.mean)[,1]) %>%
  # ungroup() %>%
  pivot_wider(names_from = type, values_from = c(ar1.sd, ar1.mean, sd.mean)) %>%
  mutate(ar1.sd.diff = scale(`ar1.sd_Fully coupled` - `ar1.sd_Mechanically decoupled`)[,1],
         ar1.mean.diff = scale(`ar1.mean_Fully coupled` - `ar1.mean_Mechanically decoupled`)[,1],
         sd.mean.diff = scale(`sd.mean_Fully coupled` - `sd.mean_Mechanically decoupled`)[,1],
         type = "Difference") %>%
  dplyr::select(lon, lat, ar1.sd.diff, ar1.mean.diff, sd.mean.diff, type) -> diff.dat


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






