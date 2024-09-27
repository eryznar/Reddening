#ERA 5 SST processing

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")


# 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
# 2) Click on "Download data" tab
# 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”

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
  nc <- nc_open("./Data/ERA5_sst_1940-2024_NEW.nc")
  
  sst <- ncvar_get(nc, "sst", verbose = F) # keep as Kelvin
  
  sst.1 <- sst[,,1,]
  sst.5 <- sst[,,2,]
  
  dim(sst.1) #241 lon, 81 lat, 770 months
  dim(sst.5) #241 lon, 81 lat, 770 months
  
  
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
  
  # Create sst data matrix
  sst.1 <- aperm(sst.1, 3:1)
  sst.5 <- aperm(sst.5, 3:1)
  
  mat_sst.1 <- t(matrix(sst.1, nrow = dim(sst.1)[1], ncol = prod(dim(sst.1)[2:3])))
  mat_sst.5 <- t(matrix(sst.5, nrow = dim(sst.5)[1], ncol = prod(dim(sst.5)[2:3])))
  
  # Convert to data frame
  data.frame(lon = lon, lat = lat,  mat_sst.1) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
    na.omit()-> sst_latlon_1960.2024 # Jan-June
  
  data.frame(lon = lon, lat = lat,  mat_sst.5) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
    na.omit()-> sst_latlon_2024 # Jul-Aug
  
  
  # Bind dataframes, calculate mean by month
  rbind(sst_latlon_1960.2024, sst_latlon_2024) %>%
    group_by(year, month, lon, lat) %>%
    reframe(mean.sst = mean(sst)) -> sst_1960.2024
  
  # Create EBS and GOA polygons, filter data into different regions
    #EBS
    xp <- cbind(ebs.x, ebs.y)
    loc=cbind(sst_1960.2024$lon, sst_1960.2024$lat)
    check <- in.poly(loc, xp=xp)
    
    cbind(check, sst_1960.2024) %>%
      filter(check == "TRUE") %>%
      dplyr::select(!check)  -> sst_EBS
    
    # Plot to check
    sst_EBS %>%
      group_by(lat, lon) %>%
      reframe(mean.sst = mean(mean.sst)) %>%
      st_as_sf(coords = c(x = "lon", y = "lat"), crs = crs.latlon) %>%
      st_transform(., map.crs) -> spat.sst.ebs
    
    cbind(st_coordinates(spat.sst.ebs), sst = spat.sst.ebs$mean.sst) -> pp
    
    layers <- akgfmaps::get_base_layers(select.region = "ebs",
                                        set.crs = map.crs)
    
    panel_extent <- data.frame(y = c(50, 70),
                               x = c(-180, -158)) %>%
      akgfmaps::transform_data_frame_crs(out.crs = map.crs)
    
    ggplot()+
      geom_point(pp, mapping = aes(X, Y, color = sst))+
      ggplot2::geom_sf(data = layers$akland,
                       fill = "grey70",
                       color = "black")+
      ggplot2::coord_sf(xlim = panel_extent$x, 
                        ylim = panel_extent$y)
    
    # Calculate mean by month
    sst_EBS %>%
      group_by(year, month) %>%
      reframe(mean.sst = mean(mean.sst)) -> sst_EBS
  
    #GOA ----
    xp <- cbind(goa.x, goa.y)
    loc=cbind(sst_1960.2024$lon, sst_1960.2024$lat)
    check <- in.poly(loc, xp=xp)
    
    cbind(check, sst_1960.2024) %>%
      filter(check == "TRUE")  %>%
      dplyr::select(!check) -> sst_GOA
    
    
    # Plot to check
    sst_GOA %>%
      group_by(lat, lon) %>%
      reframe(mean.sst = mean(mean.sst)) %>%
      st_as_sf(coords = c(x = "lon", y = "lat"), crs = crs.latlon) %>%
      st_transform(., map.crs) -> spat.sst.goa
    
    cbind(st_coordinates(spat.sst.goa), sst = spat.sst.goa$mean.sst) -> pp
    
    layers <- akgfmaps::get_base_layers(select.region = "ebs",
                                        set.crs = map.crs)
    
    panel_extent <- data.frame(y = c(38, 70),
                               x = c(-130, 125)) %>%
      akgfmaps::transform_data_frame_crs(out.crs = map.crs)
    
    ggplot()+
      geom_point(pp, mapping = aes(X, Y, color = sst))+
      ggplot2::geom_sf(data = layers$akland,
                       fill = "grey70",
                       color = "black")+
      ggplot2::coord_sf(xlim = panel_extent$x, 
                        ylim = panel_extent$y)
  
    # Calculate mean by month
    sst_GOA %>%
    group_by(year, month) %>%
      reframe(mean.sst = mean(mean.sst)) -> sst_GOA
    
    # Processing for non-stationary dynamics
    # now get monthly means weighted by area, using an arithmetic mean
    weight <- sqrt(cos(sst_GOA$lat*pi/180))
    
    sst_GOA %>%
      mutate(sst_mu = apply(mean.sst, 1, function(x) weighted.mean(x, weight, na.rm = T)))
    
    SST.mu <- apply(SST, 1, function(x) weighted.mean(x,weight, na.rm=T))
    
    # now separate out winter
    m <- months(d)
    yr <- as.numeric(as.character(years(d)))
    win <- c("Nov", "Dec", "Jan", "Feb", "Mar")
    
    mean <- data.frame(year=yr, month=m, mean=SST.mu)
    mean$win.yr <- mean$year
    mean$win.yr[mean$month %in% c("Nov", "Dec")] <- mean$win.yr[mean$month %in% c("Nov", "Dec")] + 1
    
    win.mean <- mean[mean$month %in% win,]
    
    SST.win <- tapply(win.mean$mean, win.mean$win.yr, mean)
    
    salm.cov$SST.win <- SST.win[match(salm.cov$year,(names(SST.win)))]
    
    
    
    
    
    sst_GOA %>%
      filter(month %in% c("Nov", "Dec", "Jan", "Feb", "Mar")) %>%
      mutate(year = case_when((month %in% c("Nov", "Dec")) ~ as.numeric(as.character(year)) + 1,
                              TRUE ~ as.numeric(as.character(year)))) %>% # center winter months on year of Nov-Dec
      group_by(year) %>%
      reframe(mean.sst = mean(mean.sst) - 273.15) -> winter.GOA.sst
    
    winter.GOA.sst.3 <- rollapply(winter.GOA.sst$mean.sst, 3, mean, na.rm = T, fill = NA) #3-year
    cbind(winter.GOA.sst$year, winter.GOA.sst.3) -> winter.GOA.sst
  
  
  write.csv(sst_1960.2024, "./Data/sst_1960.2024_ALL.csv")
  write.csv(sst_EBS, "./Data/sst_EBS.csv")
  write.csv(sst_GOA, "./Data/sst_GOA.csv")
  