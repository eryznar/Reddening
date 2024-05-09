lme <- st_read("./Data/LME shapefiles/LMEs66.shp")
north40 <- c(1, 2, 8, 9, 18, 19, 20, 21, 22, 23, 24, 51, 52, 53, 54, 56, 57, 58, 59, 60, 62, 63, 64)
subset(lme, lme$ARCTIC=="Arctic") -> arcticlme
subset(lme, lme$LME_NUMBER %in% north40) -> north40lme


### PROCESS SST DATA ------------------------------------------------------------------------------------------------------
nc.early <- nc_open("./Data/ERA5_sst_1960-1992.nc")
nc.late <- nc_open("./Data/ERA5_sst_1993-2024.nc")

nc <- nc.early

sst <- ncvar_get(nc, "sst", verbose = F) # keep as Kelvin
# 
# sst.1 <- sst[,,1,]
# sst.5 <- sst[,,2,]

# dim(sst.1) #1440 lon, 721 lat, 772 months
# dim(sst.5) #1440 lon, 721 lat, 772 months

dim(sst) #1440 lon, 721 lat, 396 months

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

# #Create sst data matrix
# sst.1 <- aperm(sst.1, 3:1)
# sst.5 <- aperm(sst.5, 3:1)

sst <- aperm(sst, 3:1)

mat_sst <- t(matrix(sst, nrow = dim(sst)[1], ncol = prod(dim(sst)[2:3])))
# mat_sst.5 <- t(matrix(sst.5, nrow = dim(sst.5)[1], ncol = prod(dim(sst.5)[2:3])))


data.frame(lon = lon, lat = lat,  mat_sst) %>%
  pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
  mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
  na.omit()-> sst_latlon_1960_1992

# 1993-2024 --------
nc.late <- nc_open("./Data/ERA5_sst_1993-2024.nc")

nc <- nc.late

sst <- ncvar_get(nc, "sst", verbose = F) # keep as Kelvin

sst.1 <- sst[,,1,]
sst.5 <- sst[,,2,]

dim(sst.1) #1440 lon, 721 lat, 376 months
dim(sst.5) #1440 lon, 721 lat, 376 months


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
  na.omit()-> sst_latlon_1993_2023

data.frame(lon = lon, lat = lat,  mat_sst.5) %>%
  pivot_longer(cols = c(3:ncol(.)), names_to = "month", values_to = "sst") %>%
  mutate(month = rep(m, nrow(.)/length(m)), year = rep(yr, nrow(.)/length(yr)), sst = sst) %>%
  na.omit()-> sst_latlon_2024 ## July-October


#bind with earlier sst data, subset to BB extent, calculate mean by month
rbind(sst_latlon_1960_2023, sst_latlon_2024) %>%
  group_by(year, month, lon, lat) %>%
  reframe(mean.sst = mean(sst)) %>%
  mutate(mean.sst = mean.sst - 273.15) -> sst_1960.2024

