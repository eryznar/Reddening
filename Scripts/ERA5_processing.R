#ERA 5 SST processing

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

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
  
  #Create sst data matrix
  sst.1 <- aperm(sst.1, 3:1)
  sst.5 <- aperm(sst.5, 3:1)
  
  mat_sst.1 <- t(matrix(sst.1, nrow = dim(sst.1)[1], ncol = prod(dim(sst.1)[2:3])))
  mat_sst.5 <- t(matrix(sst.5, nrow = dim(sst.5)[1], ncol = prod(dim(sst.5)[2:3])))
  
  
  data.frame(lon = lon, lat = lat,  mat_sst.1) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
    na.omit()-> sst_latlon_1960_2023 ## July-October
  
  data.frame(lon = lon, lat = lat,  mat_sst.5) %>%
    pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
    mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
    na.omit()-> sst_latlon_2024 ## July-October
  
  
  #bind with earlier sst data, subset to BB extent, calculate mean by month
  rbind(sst_latlon_1960_2023, sst_latlon_2024) %>%
    group_by(year, month, lon, lat) %>%
    reframe(mean.sst = mean(sst)) %>%
    mutate(mean.sst = mean.sst - 273.15) -> sst_1960.2024
  
  # Create EBS and GOA polygons, filter data into different regions
    #EBS
    xp <- cbind(ebs.x, ebs.y)
    loc=cbind(sst_1960.2024$lon, sst_1960.2024$lat)
    check <- in.poly(loc, xp=xp)
    
    cbind(check, sst_1960.2024) %>%
      filter(check == "TRUE") %>%
      dplyr::select(!check)  %>%
      group_by(year, month) %>%
      reframe(mean.sst = mean(mean.sst)) -> sst_EBS
  
    #GOA
    xp <- cbind(goa.x, goa.y)
    loc=cbind(sst_1960.2024$lon, sst_1960.2024$lat)
    check <- in.poly(loc, xp=xp)
    
    cbind(check, sst_1960.2024) %>%
      filter(check == "TRUE")  %>%
      dplyr::select(!check) %>%
      group_by(year, month) %>%
      reframe(mean.sst = mean(mean.sst)) -> sst_GOA
  
  
  
  write.csv(sst_1960.2024, "./Data/sst_1960.2024_ALL.csv")
  write.csv(sst_EBS, "./Data/sst_EBS.csv")
  write.csv(sst_GOA, "./Data/sst_GOA.csv")
  