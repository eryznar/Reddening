# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# several packages that are helpful for this process have been removed from CRAN, 
# apparently because dependency PCICt was itself removed

# these downloads may require package "remotes"

# download dependency PCICt
# remotes::install_version("PCICt", version = "0.5-4.3", repos = "http://cran.r-project.org")

# download ncdf4.helper
# install_version("ncdf4.helpers", version = "0.3-7")

# install.packages("https://cran.rstudio.com/src/contrib/Archive/ClimProjDiags/ClimProjDiags_0.3.3.tar.gz", 
#                  repos = NULL, type = "source")

# install s2dv
# remotes::install_version("s2dv", version = "2.1.0", repos = "https://cloud.r-project.org")

library(ncdf4)
library(ncdf4.helpers)
library(s2dv)


# Get map layers
# mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

# yrs <- 1850:2013 # these are years that are similar across SLP and SST models

# Color palettes
library(oce)
new.col <- oceColorsPalette(64)

cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

theme_set(theme_bw())

# Identify folders to download files from
fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM/"
fcm.slp.dir <- "./Data/CESM2 ensemble/SLP/FCM/"

# fcm.sst.dir <- paste0(dir, "./Data/CESM2 ensemble/SST/FCM/") #FCM SST

# 
# fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
# mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP


## Start with SST EOF/PC-----------

# Extract time info for processing below (same across files)
files <- list.files(fcm.sst.dir, full.names = TRUE)


# this is clunky - get the correct time variable for CESM using ncdf4.helper
nc.sst <- nc_open(files[1])
d <- nc.get.time.series(nc.sst)

# Split each string at "-" and put results into a matrix
date_parts <- str_split_fixed(d, "-", n = 3)

# Extract year and month as separate vectors
year <- as.numeric(date_parts[, 1])
month <- as.numeric(date_parts[, 2])

# reload the NetCDF file
nc.sst <- tidync(files[1])

# Apply spatial filter on lat and lon using hyper_filter

# and also filter for 1949-2014 (this will capture winter 1950-2014)

wanted <- which(year >= 1949 & year <= 2014)

# limit d, month, year to same years
d <- d[year >= 1949 & year <= 2014]
month <- month[year >= 1949 & year <= 2014]
year <- year[year >= 1949 & year <= 2014]

SST <- nc.sst %>%
  hyper_filter(lat = between(lat, 20, 70),
               lon = between(lon, 120, 250),
               time = index %in% wanted) %>% 
  hyper_array(select_var = "SST")



# Filter for the selected lat/lon after hyper_filter and pull their values as x,y

lon_tbl <- attr(SST, "transforms")$lon

x <- lon_tbl %>% 
  filter(selected) %>% 
  pull(lon)


lat_tbl <- attr(SST, "transforms")$lat

y <- lat_tbl %>% 
  filter(selected) %>% 
  pull(lat)

x; y

SST <- SST$SST

dim(SST) 
names(dim(SST)) <- c("lon", "lat", "time")
# 105 lon x 53 lat x 780 months


# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  


# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))


dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# and plot 
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(125,255))
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# AR1 values will be calculated on detrended monthly anomalies

# # identify columns in SST matrix corresponding to land
# land <- is.na(colMeans(SST)) 
# 
# # For analysis, we only use the columns of the matrix with non-missing values:
# X <- SST[,!land] 

# function to compute monthly means for a single time series
f <- function(x) tapply(x, month, mean, na.rm = TRUE)  

mu <- apply(SST, 2, f)	# compute monthly means for each time series (cell)

mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

anom <- SST - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
xx = seq(1,nrow(anom))

for(i in 1:ncol(anom)) {

  
  if(is.na(sum(anom[,i]))) 
    {anom.detr[,i] <- NA} 
  else
    {anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
    }
}

# limit to winter (NDJFM) 
# (peak PDO months from Newman et al. 2016)

d.NDJFM <- d[month %in% c(11,12,1:3)]
yr.NDJFM <- year[month %in% c(11,12,1:3)]
anom.detr.NDJFM <- anom.detr[month %in% c(11,12,1:3),]

# set up winter year
win.month <- month[month %in% c(11,12,1:3)]
win.yr <- if_else(win.month %in% 11:12, yr.NDJFM + 1, yr.NDJFM)

## now limit to EBS and GOA

# EBS polygon
ebs.x <- c(187, 187, 203, 203, 191) 
ebs.y <- c(53, 61, 61, 57.5, 53)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# set up SST objects for each region
#first, ebs
ebs.sst <- anom.detr.NDJFM

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# now, goa
goa.sst <- anom.detr.NDJFM

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# get annual NDJFM means for each
ebs.NDJFM.mean <- tapply(rowMeans(ebs.sst, na.rm = T), win.yr, mean)
goa.NDJFM.mean <- tapply(rowMeans(goa.sst, na.rm = T), win.yr, mean)

# limit to 1950-2014 (complete winters)
ebs.NDJFM.mean <- ebs.NDJFM.mean[names(ebs.NDJFM.mean) %in% 1950:2014]
goa.NDJFM.mean <- goa.NDJFM.mean[names(goa.NDJFM.mean) %in% 1950:2014]

# plot to check
plot.dat <- data.frame(win.year = 1950:2014,
                       EBS_anomaly = ebs.NDJFM.mean,
                       GOA_anomaly = goa.NDJFM.mean) %>%
  pivot_longer(cols = -win.year)

ggplot(plot.dat, aes(win.year, value, color = name)) + 
  geom_hline(yintercept = 0) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)])

# looks right enough for regional model output

# finally. AR(1) for 15-year rolling windows
fcm.ar <- data.frame()

for(i in 8:(length(ebs.NDJFM.mean)-7)){
  # i <- 8
  temp.ebs <- ebs.NDJFM.mean[(i-7):(i+7)]
  temp.goa <- goa.NDJFM.mean[(i-7):(i+7)]
  
  fcm.ar <- rbind(fcm.ar,
                  data.frame(winter.year = 1949+i,
                  model = "FCM",
                  ebs.ar = acf(temp.ebs, plot=FALSE)$acf[2],
                  goa.ar = acf(temp.goa, plot=FALSE)$acf[2]))
  
}

# plot to check 

plot <- fcm.ar %>%
  pivot_longer(cols = c(-winter.year, -model))

ggplot(plot, aes(winter.year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)])

### now repeat with MDM and compare!! -----

mdm.sst.dir <- "Data/CESM2 ensemble/SST/MDM/" #MDM SST

# Extract time info for processing below (same across files)
files <- list.files(mdm.sst.dir, full.names = TRUE)

# this is clunky - get the correct time variable for CESM using ncdf4.helper
nc.sst <- nc_open(files[1])
d <- nc.get.time.series(nc.sst)

# Split each string at "-" and put results into a matrix
date_parts <- str_split_fixed(d, "-", n = 3)

# Extract year and month as separate vectors
year <- as.numeric(date_parts[, 1])
month <- as.numeric(date_parts[, 2])

# reload the NetCDF file
nc.sst <- tidync(files[1])

# Apply spatial filter on lat and lon using hyper_filter

# and also filter for 1949-2014 (this will capture winter 1950-2014)

wanted <- which(year >= 1949 & year <= 2014)

# limit d, month, year to same years
d <- d[year >= 1949 & year <= 2014]
month <- month[year >= 1949 & year <= 2014]
year <- year[year >= 1949 & year <= 2014]

SST <- nc.sst %>%
  hyper_filter(lat = between(lat, 20, 70),
               lon = between(lon, 120, 250),
               time = index %in% wanted) %>% 
  hyper_array(select_var = "SST")



# Filter for the selected lat/lon after hyper_filter and pull their values as x,y

lon_tbl <- attr(SST, "transforms")$lon

x <- lon_tbl %>% 
  filter(selected) %>% 
  pull(lon)


lat_tbl <- attr(SST, "transforms")$lat

y <- lat_tbl %>% 
  filter(selected) %>% 
  pull(lat)

x; y

SST <- SST$SST

dim(SST) 
names(dim(SST)) <- c("lon", "lat", "time")
# 105 lon x 53 lat x 780 months


# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  


# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))


dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# and plot 
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(125,255))
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# AR1 values will be calculated on detrended monthly anomalies

# # identify columns in SST matrix corresponding to land
# land <- is.na(colMeans(SST)) 
# 
# # For analysis, we only use the columns of the matrix with non-missing values:
# X <- SST[,!land] 

# function to compute monthly means for a single time series
f <- function(x) tapply(x, month, mean, na.rm = TRUE)  

mu <- apply(SST, 2, f)	# compute monthly means for each time series (cell)

mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

anom <- SST - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
xx = seq(1,nrow(anom))

for(i in 1:ncol(anom)) {
  
  
  if(is.na(sum(anom[,i]))) 
  {anom.detr[,i] <- NA} 
  else
  {anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
  }
}

# limit to winter (NDJFM) 
# (peak PDO months from Newman et al. 2016)

d.NDJFM <- d[month %in% c(11,12,1:3)]
yr.NDJFM <- year[month %in% c(11,12,1:3)]
anom.detr.NDJFM <- anom.detr[month %in% c(11,12,1:3),]

# set up winter year
win.month <- month[month %in% c(11,12,1:3)]
win.yr <- if_else(win.month %in% 11:12, yr.NDJFM + 1, yr.NDJFM)

## now limit to EBS and GOA

# EBS polygon
ebs.x <- c(187, 187, 203, 203, 191) 
ebs.y <- c(53, 61, 61, 57.5, 53)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

# set up SST objects for each region
#first, ebs
ebs.sst <- anom.detr.NDJFM

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# now, goa
goa.sst <- anom.detr.NDJFM

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# get annual NDJFM means for each
ebs.NDJFM.mean <- tapply(rowMeans(ebs.sst, na.rm = T), win.yr, mean)
goa.NDJFM.mean <- tapply(rowMeans(goa.sst, na.rm = T), win.yr, mean)

# limit to 1950-2014 (complete winters)
ebs.NDJFM.mean <- ebs.NDJFM.mean[names(ebs.NDJFM.mean) %in% 1950:2014]
goa.NDJFM.mean <- goa.NDJFM.mean[names(goa.NDJFM.mean) %in% 1950:2014]

# plot to check
plot.dat <- data.frame(win.year = 1950:2014,
                       EBS_anomaly = ebs.NDJFM.mean,
                       GOA_anomaly = goa.NDJFM.mean) %>%
  pivot_longer(cols = -win.year)

ggplot(plot.dat, aes(win.year, value, color = name)) + 
  geom_hline(yintercept = 0) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)])

# looks right enough for regional model output

# finally. AR(1) for 15-year rolling windows
mdm.ar <- data.frame()

for(i in 8:(length(ebs.NDJFM.mean)-7)){
  # i <- 8
  temp.ebs <- ebs.NDJFM.mean[(i-7):(i+7)]
  temp.goa <- goa.NDJFM.mean[(i-7):(i+7)]
  
  mdm.ar <- rbind(mdm.ar,
                  data.frame(winter.year = 1949+i,
                             model = "MDM",
                             ebs.ar = acf(temp.ebs, plot=FALSE)$acf[2],
                             goa.ar = acf(temp.goa, plot=FALSE)$acf[2]))
  
}

# plot to check 

plot <- mdm.ar %>%
  pivot_longer(cols = c(-winter.year, -model))

ggplot(plot, aes(winter.year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)])

# plot both to compare - FCM vs MDM
plot.both <- rbind(fcm.ar, mdm.ar) %>%
  pivot_longer(cols = c(-winter.year, -model))


ggplot(plot.both, aes(winter.year, value, color = model)) +
  facet_wrap(~name) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)])

### compare with observed AR(1) ----

# load ERSST montly anomalies
ebs.obs <- read.csv("./Output/SST.winter.anom.ebs.detr.csv", row.names = 1) %>%
  filter(Year %in% 1950:2014) %>%
  mutate(winter.year = Year)

goa.obs <- read.csv("./Output/SST.winter.anom.goa.detr.csv", row.names = 1) %>%
  filter(Year %in% 1950:2014) %>%
  mutate(winter.year = Year)

# get AR(1) on 15-year moving windows for both

obs.ar <- data.frame()

for(i in 8:(nrow(ebs.obs)-7)){
   # i <- 8
  temp.ebs <- ebs.obs[(i-7):(i+7),2]
  temp.goa <- goa.obs[(i-7):(i+7),2]
  
  obs.ar <- rbind(obs.ar,
                  data.frame(winter.year = 1949+i,
                             model = "Obs",
                             ebs.ar = acf(temp.ebs, plot=FALSE)$acf[2],
                             goa.ar = acf(temp.goa, plot=FALSE)$acf[2]))
  
}

obs.ar

# plot  to compare - Obs vs FCM vs MDM
plot.all <- rbind(fcm.ar, mdm.ar, obs.ar) %>%
  pivot_longer(cols = c(-winter.year, -model))


ggplot(plot.all, aes(winter.year, value, color = model)) +
  facet_wrap(~name) +
  geom_line() +
  scale_color_manual(values = cb[c(2,4,1)])

# check for correlation between obs and FCM/MDM
ebs.cor <- plot.all %>%
  filter(name == "ebs.ar") %>%
  pivot_wider(names_from = model, values_from = value) 

cor(ebs.cor[,3:5]) # FCM a very weak correlation, MDM nothing

goa.cor <- plot.all %>%
  filter(name == "goa.ar") %>%
  pivot_wider(names_from = model, values_from = value) 

cor(goa.cor[,3:5]) # FCM a very weak correlation, MDM negatively correlated

# so that might be one indication that FCM captures the dynamics we're interested in 
# better than MDM, indicating a role of AL variability in SST reddening

# how about the SLP PC1 SD (explanatory variable) vs SST AR1 (response variable) relationship?
# does FCM capture this better than MDM?


