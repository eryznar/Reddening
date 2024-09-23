#ERA 5 SLP processing

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")

# 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
# 2) Click on "Download data" tab
# 3) Click on the product type you’d like. This script processes 

### SET STUDY REGIONS ------------------------------------------------------------------------------------------------------
# EBS:
ebs.x <- c(183, 183, 203, 203, 191)
ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)

ebs.y <- c(53, 65, 65, 57.5, 53)

# GOA: 
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.x <- ifelse(goa.x > 180, goa.x-360, goa.x)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)


### PROCESS SST DATA ------------------------------------------------------------------------------------------------------
  nc <- nc_open("./Data/ERA5_SLP_1960-2024.nc")
  
  slp <- ncvar_get(nc, "msl", verbose = F) # keep as Kelvin
  
  slp.1 <- slp[,,1,]
  slp.5 <- slp[,,2,]
  
  dim(slp.1) #161 lon, 201 lat, 776 months
  dim(slp.5) #161 lon, 201 lat, 776 months
  
  #Process
  h <- (ncvar_get(nc, "time")/24)
  d <- dates(h, origin = c(1, 1, 1900))  
  m <- months(d)
  yr <- chron::years(d)
  
  x <- ncvar_get(nc, "longitude")
  y <- ncvar_get(nc, "latitude")
  
  # Keep track of corresponding latitudes and longitudes of each column:
  lat <- rep(y, length(x))
  lon <- rep(x, each = length(y))

  # Create slp data matrix
  slp.1 <- aperm(slp.1, 3:1)
  slp.5 <- aperm(slp.5, 3:1)
  
  mat_slp.1 <- t(matrix(slp.1, nrow = dim(slp.1)[1], ncol = prod(dim(slp.1)[2:3])))
  mat_slp.5 <- t(matrix(slp.5, nrow = dim(slp.5)[1], ncol = prod(dim(slp.5)[2:3])))
  
  # Convert to data frame
  data.frame(lon = lon, lat = lat,  mat_slp.1) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "slp") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), slp = slp) %>%
    na.omit()-> slp_latlon_1960.2024 # Jan-June
  
  data.frame(lon = lon, lat = lat,  mat_slp.5) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "slp") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), slp = slp) %>%
    na.omit()-> slp_latlon_2024 # Jul-Aug
  
  
  # # Bind dataframes, calculate mean by month
  # rbind(slp_latlon_1960.2024, slp_latlon_2024) %>%
  #   group_by(year, month, lon, lat) %>%
  #   reframe(mean.slp = mean(slp)) -> slp_1960.2024
  
  # Compute monthly means for each times series (cell)
  rbind(slp_latlon_1960.2024, slp_latlon_2024) %>%
    group_by(month, lon, lat) %>%
    reframe(mu = mean(slp)) -> month.means
  
  # Join back with year data
  right_join(rbind(slp_latlon_1960.2024, slp_latlon_2024), month.means) %>% 
    na.omit() -> slp.means
  
  # Calculate anomalies
  slp.means %>%
    mutate(slp.anom = slp - mu) -> slp.means
  
  # Create EBS and GOA polygons, filter data into different regions
  #EBS
  xp <- cbind(ebs.x, ebs.y)
  loc=cbind(slp.means$lon, slp.means$lat)
  check <- in.poly(loc, xp=xp)
  
  cbind(check, slp.means) %>%
    filter(check == "TRUE", month %in% c("Jan", "Feb", "Mar", "Apr")) %>%
    dplyr::select(!check)  -> winter.slp.EBS
    # group_by(year, month) %>%
    # reframe(mean.slp = mean(mean.slp)) 
  
  # Save winter slp
  write.csv(winter.slp.EBS, "./Output/winter.slp.EBS.csv")
  
  # Get weights
  winter.slp.EBS %>%
    dplyr::select(lat) %>%
    pull() -> EBS.lat
  
  EBS.weights <- sqrt(cos(EBS.lat*pi/180))
  
  write.csv(EBS.weights, "./Output/EBS.weights.csv")
  
  #GOA
  xp <- cbind(goa.x, goa.y)
  loc=cbind(slp.means$lon, slp.means$lat)
  check <- in.poly(loc, xp=xp)
  
  cbind(check, slp.means) %>%
    filter(check == "TRUE", month %in% c("Jan", "Feb", "Mar", "Apr")) %>%
    dplyr::select(!check)  -> winter.slp.GOA
  # group_by(year, month) %>%
  # reframe(mean.slp = mean(mean.slp)) 
  
  # Save winter slp
  write.csv(winter.slp.GOA, "./Output/winter.slp.GOA.csv")
  
  # Get weights
  winter.slp.GOA %>%
    dplyr::select(lat) %>%
    pull() -> GOA.lat
  
  GOA.weights <- sqrt(cos(GOA.lat*pi/180))
  
  write.csv(GOA.weights, "./Output/GOA.weights.csv")
  
  
# CALCULATE PC1 AND RUN EOF ---------------------------------------------------------------------------
  # EBS ----
  # load slp and cell weights
  ebs.slp <- read.csv("./Output/winter.slp.EBS.csv", row.names = 1)
  
  weights <- read.csv("./Output/EBS.weights.csv", row.names = 1)
  
  # Make data wider for pca
  ebs.slp %>%
    mutate(coords = paste0("N", lat, "W", lon*-1),
           date = paste0(month, "-", year)) %>%
    dplyr::select("slp.anom", "coords", "date") %>%
    pivot_wider(values_from = slp.anom, names_from = c("coords")) %>%
    as.data.frame(.) -> ebs.slp.wide
  
  # Extract year and make date as rownames
  yr <- str_split(ebs.slp.wide$date, "-", simplify = T)[,2]
  
  rownames(ebs.slp.wide) <- ebs.slp.wide$date
  
  ebs.slp.wide %>% dplyr::select(!date) -> ebs.slp.wide
  
  # Calculate weights
  temp1 <- str_split(colnames(ebs.slp.wide), "W", simplify = T)[,1]
  X.lats <- as.numeric(str_split(temp1, "N", simplify = T)[,2])
  
  ebs.weights <- sqrt(cos(X.lats*pi/180))
    
  # Run PCA
  pca <- FactoMineR::svd.triplet(cov(ebs.slp.wide), col.w=ebs.weights) #weighting the columns
  
  pc1_slp.ebs <- as.matrix(ebs.slp.wide) %*% pca$U[,1]
  
  # and scale!
  pc1_slp.ebs <- as.vector(scale(pc1_slp.ebs))
  
  
  # and annual FMA means for this value
  pc1_slp.ebs <- tapply(pc1_slp.ebs, yr, mean)
  
  # plot to check
  plot(1960:2024, pc1_slp, type="l") # looks right re. 1976/77
  
