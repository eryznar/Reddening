# PURPOSE: to download CESM2 fcm and mdm SST/SLP model outputs and calculate ar1/SD

# Note: 

# Author: Emily Ryznar

# 
source("./Scripts/load.libs.functions.R")

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

# Extract time info for processing below
time_units <- ncmeta::nc_atts(files[1], "time") %>%
  filter(name == "units") %>%
  pull(value)

unit_parts <- str_split(time_units, " since ")[[1]]
time_unit <- unit_parts[1]
origin_date <- ymd_hms(unit_parts[2])

# Get map layers
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))


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
      group_by(lon, lat, member) %>%
      mutate(mean.cell.SST = mean(SST)) %>% # compute grid cell mean across years
      ungroup() %>%
      mutate(SSTa = SST - mean.cell.SST) %>% # compute anomalies
      group_by(lon, lat, year, member) %>%
      reframe(mean.SSTa = mean(SSTa), # calculate mean
              mean.SST = mean(SST))-> out
    
    # Detrend data and extract residuals
    resid.SST <- lm(mean.SST ~ year, out)$residuals
    resid.SSTa <- lm(mean.SSTa ~ year, out)$residuals
    
    # Bind with output df
    detrend.dat <- cbind(out, resid.SST, resid.SSTa)
    
    # Calculate grid cell AR1
    detrend.dat %>%
      group_by(lon, lat, member) %>%
      reframe(ar1.SSTa = sapply(rollapply(resid.SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
              ar1.sd = sd(ar1.SSTa),
              mean.residSSTa = mean(resid.SSTa),
              mean.SSTa = mean(mean.SSTa),
              mean.SST = mean(mean.SST)) %>%
      mutate(lon = as.numeric(lon),
             lat = as.numeric(lat)) %>%
      dplyr::select(!ar1.SSTa) %>%
      distinct() -> ar1.dat
    
    # stack processed files
    fcm.sst <- bind_rows(fcm.sst, ar1.dat)
    
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
      group_by(lon, lat, member) %>%
      mutate(mean.cell.SST = mean(SST)) %>% # compute grid cell mean across years
      ungroup() %>%
      mutate(SSTa = SST - mean.cell.SST) %>% # compute anomalies
      group_by(lon, lat, year, member) %>%
      reframe(mean.SSTa = mean(SSTa), # calculate mean
              mean.SST = mean(SST))-> out
    
    # Detrend data and extract residuals
    resid.SST <- lm(mean.SST ~ year, out)$residuals
    resid.SSTa <- lm(mean.SSTa ~ year, out)$residuals
    
    # Bind with output df
    detrend.dat <- cbind(out, resid.SST, resid.SSTa)
    
    # Calculate grid cell AR1
    detrend.dat %>%
      group_by(lon, lat, member) %>%
      reframe(ar1.SSTa = sapply(rollapply(resid.SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
              ar1.sd = sd(ar1.SSTa),
              mean.residSSTa = mean(resid.SSTa),
              mean.SSTa = mean(mean.SSTa),
              mean.SST = mean(mean.SST)) %>%
      mutate(lon = as.numeric(lon),
             lat = as.numeric(lat)) %>%
      dplyr::select(!ar1.SSTa) %>%
      distinct() -> ar1.dat
    
    
    # stack processed files
    mdm.sst <- bind_rows(mdm.sst, ar1.dat)
    
  }
  
  # Save file
  saveRDS(mdm.sst, paste0(dir, "Output/processed.mdm.sst.rda"))
  
  
  
# 3) CALCULATE/PLOT CELL-WISE AR1 SD, SSTa SD, and MEAN AR1 ACROSS ENSEMBLE ------------------------
# Calculate SD in AR1 by grid cell across ensemble members
fcm.sst <- readRDS(paste0(dir, "Output/processed.fcm.sst.rda"))

fcm.sst %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          SST.sd = sd(mean.SSTa)) %>%
  mutate(type = "FCM") -> fcm.sst.ar1sd

mdm.sst <- readRDS(paste0(dir, "Output/processed.mdm.sst.rda"))

mdm.sst %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = mean(ar1.sd),
          SST.sd = sd(mean.SSTa)) %>%
  mutate(type = "MDM") -> mdm.sst.ar1sd

# Join
plot.dat <- rbind(fcm.sst.ar1sd, mdm.sst.ar1sd)


# Plot
ggplot()+
  geom_tile(plot.dat, mapping= aes(lon, lat, fill = ar1.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)+
  facet_wrap(~type, nrow = 2)+
  #ggtitle("Ensemble SST AR1 SD")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  median(plot.dat$ar1.sd))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10))  -> ar1.sd.plot

ggplot()+
  geom_tile(plot.dat, mapping= aes(lon, lat, fill = SST.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)+
  #ggtitle("Ensemble SSTa SD")+
  facet_wrap(~type, nrow = 2)+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  median(plot.dat$SST.sd))+
  theme_bw()+
  theme(plot.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 10)) -> SSTa.sd.plot


ggsave(plot = ar1.sd.plot, "./Figures/CESM2_AR1_SD.png", width = 5, height = 5, units = "in")
ggsave(plot = SSTa.sd.plot, "./Figures/CESM2_SSTa_SD.png", width = 5, height = 5, units = "in")

