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
    }, by = .(lat, lon)]
    
    
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
    }, by = .(lat, lon)]
    
    
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
      group_by(lon, lat, win.year, member) %>% 
      reframe(mean.SLPa = mean(SLPa)) -> out # calculate mean annual SLPa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SLPa := { # output column name
      fit <- lm(mean.SLPa ~ win.year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon)]

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
      filter(month %in% c(11:12, 1:3)) %>% #filter by winter months
      group_by(lat, lon, month, member) %>%
      mutate(mean.month.SLP = mean(PSL)) %>% # compute monthly mean by grid cell and member
      ungroup() %>%
      mutate(SLPa = PSL - mean.month.SLP, # compute anomalies
             win.year = case_when((month %in% 11:12) ~ year + 1, # calculate winter year
                                  TRUE ~ year)) %>% 
      group_by(lon, lat, win.year, member) %>% 
      reframe(mean.SLPa = mean(SLPa)) -> out # calculate mean annual SLPa by grid cell
    
    
    # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
    setDT(out) # convert to data.table
    
    out[, detrended_SLPa := { # output column name
      fit <- lm(mean.SLPa ~ win.year)  # Fit linear model for each lat/lon group
      residuals(fit)           # Extract residuals as detrended values
    }, by = .(lat, lon)]
    
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
    group_by(lon, lat, member) %>%
    reframe(ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
            ar1.cv = cv(ar1.SSTa), # sd of AR1 of rolling window timeseries
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
    group_by(lon, lat, member) %>%
    reframe(ar1.sd = sd(ar1.SSTa), # sd of AR1 of rolling window timeseries
            ar1.cv = cv(ar1.SSTa), # sd of AR1 of rolling window timeseries
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
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Fully coupled") -> fcm.ar1.sd.dat2


mdm.ar1.sd.dat %>%
  filter(lat >=20 & lat<=68) %>% # isolate extra tropical N. pacific
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "Mechanically decoupled") -> mdm.ar1.sd.dat2

# Join
plot.dat <- rbind(fcm.ar1.sd.dat2, mdm.ar1.sd.dat2)

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







# 4) EVALUATE PERIODS OF HIGH SLP VARIABILITY ----
slp.fcm <- readRDS(paste0(dir, "Output/FCM_winterSLPa_ar1sd.rda"))

mems <- unique(slp.fcm$member)


slp.fcm %>%
  filter(member == mems[16]) %>% # isolating one ensemble member
  group_by(win.year) %>%
  reframe(sd = mean(sd.SLPa)) %>%
  na.omit()-> tt

# Seems like 1984-1921 has high AL variability and 1922-1947 has low
ggplot(tt %>% filter(win.year > 1857), aes(win.year, sd))+
  geom_point()+
  geom_line() + 
  scale_x_continuous(breaks = seq(min(tt$win.year), max(tt$win.year), by = 10))+
  theme_bw() +
  ylab("SLP standard deviation")+
  ggtitle("AL high activity area winter SLPa")

# # By ensemble member
# slp.fcm %>%
#   group_by(win.year, member) %>%
#   reframe(sd = mean(sd.SLPa)) %>%
#   na.omit()-> tt
# 
# 
# ggplot(tt %>% filter(win.year > 1857, member %in% mems[1:10]), aes(win.year, sd, group = member))+
#   #geom_point()+
#   geom_line() + 
#   scale_x_continuous(breaks = seq(min(tt$win.year), max(tt$win.year), by = 10))+
#   theme_bw() +
#   ylab("SLP standard deviation")+
#   ggtitle("AL high activity area winter SLPa")+
#   theme(legend.position = "none")

highyrs <- 1884:1918
lowyrs <- 1919:1947

# Process FCM and MDM SST in winter for comparison in high and low SLP var periods ----
files <- list.files(fcm.sst.dir, full.names = TRUE)

time_units <- ncmeta::nc_atts(files[1], "time") %>% 
  filter(name == "units") %>%
  pull(value)

unit_parts <- str_split(time_units, " since ")[[1]]
time_unit <- unit_parts[1]
origin_date <- ymd_hms(unit_parts[2])

# FCM SST  
files <- list.files(fcm.sst.dir, full.names = TRUE)

fcm.win.sst <- tibble()

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
           member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
    group_by(lat, lon, month, member) %>%
    mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
    ungroup() %>%
    mutate(SSTa = SST - mean.month.SST,
           win.year = case_when((month %in% 11:12) ~ year + 1,
                                TRUE ~ year)) %>% # compute anomalies
    filter(win.year %in% c(highyrs, lowyrs),
           month %in% c(11:12, 1:3)) %>%
    mutate(period = case_when((win.year %in% highyrs) ~ "high",
                              TRUE ~ "low")) %>%
    group_by(lon, lat, win.year, member, period) %>% 
    reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell
  
  
  # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
  setDT(out) # convert to data.table
  
  out[, detrended_winSSTa := { # output data table column
    fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
    residuals(fit)           # Extract residuals as detrended values
  }, by = .(lat, lon, period)]
  
  
  # Calculate rolling window AR1 within data.table
  out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
    acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
    return(acf_result$acf[2])
  }, align = "center"), by = c("lon", "lat", "member", "period")]
  
  # Calculate rolling window SD within data.table
  out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd, 
                              align = "center"), by = c("lon", "lat", "member", "period")]
  
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
           member = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
    group_by(lat, lon, month, member) %>%
    mutate(mean.month.SST = mean(SST)) %>% # compute monthly mean by grid cell and member
    ungroup() %>%
    mutate(SSTa = SST - mean.month.SST,
           win.year = case_when((month %in% 11:12) ~ year + 1,
                                TRUE ~ year)) %>% # compute anomalies
    filter(win.year %in% c(highyrs, lowyrs),
           month %in% c(11:12, 1:3)) %>%
    mutate(period = case_when((win.year %in% highyrs) ~ "high",
                              TRUE ~ "low")) %>%
    group_by(lon, lat, win.year, member, period) %>% 
    reframe(mean.winSSTa = mean(SSTa))-> out # calculate mean annual winter SSTa by grid cell
  
  
  # Detrend data and extract residuals, using data.table package (AWESOME!) to speed things up
  setDT(out) # convert to data.table
  
  out[, detrended_winSSTa := { # output data table column
    fit <- lm(mean.winSSTa ~ win.year)  # Fit linear model for each lat/lon group and period(?)
    residuals(fit)           # Extract residuals as detrended values
  }, by = .(lat, lon, period)]
  
  
  # Calculate rolling window AR1 within data.table
  out[, ar1.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = function(x) {
    acf_result <- acf(x, lag.max = 1, plot = FALSE, na.action = na.pass)
    return(acf_result$acf[2])
  }, align = "center"), by = c("lon", "lat", "member", "period")]
  
  # Calculate rolling window SD within data.table
  out[, sd.winSSTa := frollapply(detrended_winSSTa, n = 15, FUN = sd, 
                                 align = "center"), by = c("lon", "lat", "member", "period")]
  
  # stack processed files
  mdm.win.sst <- bind_rows(mdm.win.sst, out)
  
}

setDT(mdm.win.sst)

# Save file
saveRDS(mdm.win.sst, paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda"))


# Calculate and plot cell-wise AR1, SD, mean AR1 ----
fcm.sst <- readRDS(paste0(dir, "Output/FCM_winterSSTa_ar1sd.rda"))

# Calculate grid cell AR1
fcm.sst %>%
  na.omit() %>%
  group_by(lon, lat, member, period) %>%
  reframe(ar1.sd = sd(ar1.winSSTa), # sd of AR1 of rolling window timeseries
          ar1.cv = (sd(ar1.winSSTa)/mean(ar1.winSSTa))*100, # sd of AR1 of rolling window timeseries
          ar1.mean = mean(ar1.winSSTa), # mean AR1 of rolling window timeseries
          sd.mean = mean(sd.winSSTa)) %>% # mean SD of rolling window timeseries
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat)) %>%
  distinct() -> fcm.ar1.sd.dat

# MDM SST ----
mdm.sst <- readRDS(paste0(dir, "Output/MDM_winterSSTa_ar1sd.rda"))

# Calculate grid cell AR1
mdm.sst %>%
  na.omit() %>%
  group_by(lon, lat, member, period) %>%
  reframe(ar1.sd = sd(ar1.winSSTa), # sd of AR1 of rolling window timeseries
          ar1.cv = (sd(ar1.winSSTa)/mean(ar1.winSSTa))*100, # sd of AR1 of rolling window timeseries
          ar1.mean = mean(ar1.winSSTa), # mean AR1 of rolling window timeseries
          sd.mean = mean(sd.winSSTa)) %>% # mean SD of rolling window timeseries
  mutate(lon = as.numeric(lon),
         lat = as.numeric(lat)) %>%
  distinct() -> mdm.ar1.sd.dat



# Calculate mean SD AR1 and AR1 by grid cell across ensemble members
fcm.ar1.sd.dat %>%
  group_by(lon, lat, period) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.cv = mean(ar1.cv),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "FCM") -> fcm.ar1.sd.dat2


mdm.ar1.sd.dat %>%
  group_by(lon, lat, period) %>%
  reframe(ar1.sd = mean(ar1.sd),
          ar1.cv = mean(ar1.cv),
          ar1.mean = mean(ar1.mean),
          sd.mean = mean(sd.mean)) %>%
  mutate(type = "MDM") -> mdm.ar1.sd.dat2

# Join
plot.dat <- rbind(fcm.ar1.sd.dat2, mdm.ar1.sd.dat2)

# Calculate differences
plot.dat %>%
  group_by(period, lon, lat) %>%
  reframe(ar1.sd.diff = ar1.sd[type == `Fully coupled`] - ar1.sd[type == `Mechanically decoupled`]) ,
          ar1.cv.diff = scale(ar1.cv[type == "Fully coupled"] - ar1.cv[type == "Mechanically decoupled"])[,1],
          ar1.mean.diff = scale(ar1.mean[type == "Fully coupled"] - ar1.mean[type == "Mechanically decoupled"])[,1],
          sd.mean.diff = scale(sd.mean[type == "Fully coupled"] - sd.mean[type == "Mechanically decoupled"])[,1]) %>%
  mutate(type = "Difference") -> diff.dat

plot.dat %>%
  pivot_wider(names_from = c("type", "period"), values_from = c("ar1.sd", "ar1.cv", "ar1.mean", "sd.mean")) %>%
  mutate(ar1.sd.high = scale(ar1.sd_FCM_high - ar1.sd_MDM_high)[,1],
         ar1.sd.low = scale(ar1.sd_FCM_low - ar1.sd_MDM_low)[,1],
         ar1.cv.high = scale(ar1.cv_FCM_high - ar1.cv_MDM_high)[,1],
         ar1.cv.low = scale(ar1.cv_FCM_low - ar1.cv_MDM_low)[,1],
         ar1.mean.high = scale(ar1.mean_FCM_high - ar1.mean_MDM_high)[,1],
         ar1.mean.low = scale(ar1.mean_FCM_low - ar1.mean_MDM_low)[,1],
         sd.mean.high = scale(sd.mean_FCM_high - sd.mean_MDM_high)[,1],
         sd.mean.low = scale(sd.mean_FCM_low - sd.mean_MDM_low)[,1]) %>%
  dplyr::select(lon, lat, ar1.sd.high, ar1.sd.low, ar1.cv.high, ar1.cv.low, ar1.mean.high, ar1.mean.low, sd.mean.high,
                sd.mean.low) -> diff.dat

diff.dat %>%
  pivot_longer(!c(lon, lat)) %>%
  mutate(name = rep(c("ar1.sd", "ar1.sd", "ar1.cv", "ar1.cv", "ar1.mean", "ar1.mean", "sd.mean", "sd.mean"),
                     nrow(.)/8),
         period = rep(c("high", "low"), nrow(.)/2)) %>%
  mutate(type = "Difference") -> diff.dat2

# Plot AR1 SD and difference
ggplot()+
  geom_tile(plot.dat %>% filter(period == "high"), mapping= aes(lon, lat, fill = ar1.sd))+
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

dat <- diff.dat2 %>% filter(name == "ar1.sd", period == "high")
ggplot()+
  geom_tile(dat, mapping= aes(lon, lat, fill = value))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 68), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type)+
  xlab("Latitude")+
  ylab("Longitude")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  0,
                       name = "FCM-MDM diff",
                       limits = c(max(abs(dat$value))*-1, max(abs(dat$value))))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 8),
        axis.title = element_text(size = 10)) -> ar1.sd.diffplot


ar1.sd.plot + ar1.sd.diffplot + plot_layout(nrow = 2, ncol = 1, byrow = TRUE, 
                                            widths = c(1, 1), heights = c(1, 0.5),
                                            axes = "collect")

#ggsave("./Figures/CESM2_AR1_SD.png", height= 7, width = 5, units = "in")

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






