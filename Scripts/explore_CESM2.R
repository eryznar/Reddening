source("./Scripts/load.libs.functions.R")

# DOWNLOAD CESM2 MODEL OUTPUTS ----------------------------------------------
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

# LOAD AND STACK CESM MODEL OUTPUTS --------------------------------
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
           model = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
    group_by(lon, lat, model) %>%
    mutate(mean.cell.SST = mean(SST)) %>% # compute grid cell mean across years
    ungroup() %>%
    mutate(SSTa = SST - mean.cell.SST) %>% # compute anomalies
    group_by(lon, lat, year, model) %>%
    reframe(mean.SSTa = mean(SSTa), # calculate mean
            mean.SST = mean(SST))-> out
  
  # Detrend data and extract residuals
  resid.SST <- lm(mean.SST ~ year, out)$residuals
  resid.SSTa <- lm(mean.SSTa ~ year, out)$residuals
  
  # Bind with output df
  detrend.dat <- cbind(out, resid.SST, resid.SSTa)
  
  # Calculate grid cell AR1
  detrend.dat %>%
    group_by(lon, lat, model) %>%
    reframe(ar1.SSTa = sapply(rollapply(resid.SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
            mean.residSSTa = mean(resid.SSTa),
            mean.SSTa = mean(mean.SSTa),
            mean.SST = mean(mean.SST)) %>%
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat)) %>%
    rename(member = model) -> ar1.dat
  
  # stack processed files
  fcm.sst <- bind_rows(fcm.sst, ar1.dat)
  
}
rbind(fcm.sst, fcm.sst2) -> fcm.sst.new

# Calculate SD in AR1 by grid cell across ensemble members
fcm.sst.new %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = sd(ar1.SSTa),
          SST.sd = sd(mean.SSTa),
          mean.ar1 = mean(ar1.SSTa)) -> fcm.sst.ar1sd

# Plot
ggplot()+
  geom_tile(fcm.sst.ar1sd, mapping= aes(lon, lat, fill = ar1.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)+
  ggtitle("Ensemble AR1 SD")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  median(fcm.sst.ar1sd$ar1.sd))

ggplot()+
  geom_tile(fcm.sst.ar1sd, mapping= aes(lon, lat, fill = SST.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)+
  ggtitle("Ensemble SSTa SD")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  median(fcm.sst.ar1sd$SST.sd))

ggplot()+
  geom_tile(fcm.sst.ar1sd, mapping= aes(lon, lat, fill = mean.ar1))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)+
  ggtitle("Ensemble mean AR1")+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white", 
                       midpoint=  median(fcm.sst.ar1sd$mean.ar1))

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
           model = substr(files[ii], 81, 88)) %>% # extracting ensemble member #
    group_by(lon, lat, model) %>%
    mutate(mean.cell.SST = mean(SST)) %>% # compute grid cell mean across years
    ungroup() %>%
    mutate(SSTa = SST - mean.cell.SST) %>% # compute anomalies
    group_by(lon, lat, year, model) %>%
    reframe(mean.SSTa = mean(SSTa), # calculate mean
            mean.SST = mean(SST))-> out
  
  # Detrend data and extract residuals
  resid.SST <- lm(mean.SST ~ year, out)$residuals
  resid.SSTa <- lm(mean.SSTa ~ year, out)$residuals
  
  # Bind with output df
  detrend.dat <- cbind(out, resid.SST, resid.SSTa)
  
  # Calculate grid cell AR1
  detrend.dat %>%
    group_by(lon, lat, model) %>%
    reframe(ar1.SSTa = sapply(rollapply(resid.SSTa, width = 15, FUN = acf, lag.max = 1, plot = FALSE)[,1], "[[",2),
            mean.residSSTa = mean(resid.SSTa),
            mean.SSTa = mean(mean.SSTa),
            mean.SST = mean(mean.SST)) %>%
    mutate(lon = as.numeric(lon),
           lat = as.numeric(lat)) %>%
    rename(member = model) -> ar1.dat
  
  # stack processed files
  mdm.sst <- bind_rows(mdm.sst, ar1.dat)
  
}

# Calculate SD in AR1 by grid cell across ensemble members
mdm.sst %>%
  group_by(lon, lat) %>%
  reframe(ar1.sd = sd(ar1)) -> mdm.sst.ar1sd

# Plot
ggplot()+
  geom_tile(mdm.sst.ar1sd, mapping= aes(lon, lat, fill = ar1.sd))+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(0, 75), xlim = c(125, 255), expand = FALSE)

# Calculate SSTa SD by grid cell across ensemble members


### PROCESS SLPa for NORTH PACIFIC REGION -----------------------------------------
# first, load data
nc.slp <- nc_open(paste0(dir, "Data/hawaii_soest_f19d_3925_d70b_1322_e90d_09e0NEW.nc"))

# process SLP data - first, extract dates
raw <- ncvar_get(nc.slp, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract coordinates
x <- ncvar_get(nc.slp, "longitude")
y <- ncvar_get(nc.slp, "latitude")

# extract data
SLP <- ncvar_get(nc.slp, "slp", verbose = F)

# Change data to a matrix
SLP <- aperm(SLP, 3:1)  

# Change to matrix
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Get lat/long vectors and add names to SLP matrix
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# Filter to area of interest
poly.x <- c(130, 130, 250, 250, 130) 
poly.y <- c(5, 70, 70, 5, 5)