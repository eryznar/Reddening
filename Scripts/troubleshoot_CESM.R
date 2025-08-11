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
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

yrs <- 1850:2013 # these are years that are similar across SLP and SST models

# Color palette
library(oce)

new.col <- oceColorsPalette(64)

theme_set(theme_bw())

# Identify folders to download files from
fcm.sst.dir <- "./Data/CESM2 ensemble/SST/FCM/"
fcm.slp.dir <- "./Data/CESM2 ensemble/SLP/FCM/"

# fcm.sst.dir <- paste0(dir, "./Data/CESM2 ensemble/SST/FCM/") #FCM SST
# mdm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/MDM/") #MDM SST
# 
# fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
# mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP


## Start with SST EOF/PC-----------

# Extract time info for processing below (same across files)
files <- list.files(fcm.sst.dir, full.names = TRUE)

# Calculate EOF on SST data
# first, load data
nc.sst <- nc_open(files[1]) # Isolate one CESM member

# d <- nc.sst



# this is clunky - get the correct time variable for CESM using ncdf4.helper
nc.sst <- nc_open(files[1])
d <- nc.get.time.series(nc.sst)
# 
# names(dim(SST)) <- c("lon", "lat", "time")

# now trim SST to 1950:2014 (when we expect EOF to be the PDO!)
# Split each string at "-" and put results into a matrix
date_parts <- str_split_fixed(d, "-", n = 3)

# Extract year and month as separate vectors
year <- as.numeric(date_parts[, 1])
month <- as.numeric(date_parts[, 2])

# reload the NetCDF file
nc.sst <- tidync(files[1])

# Apply spatial filter on lat and lon using hyper_filter

# and also filter for 1950-2014

wanted <- which(year >= 1950 & year <= 2014)

# limit d, month, year to same years
d <- d[year >= 1950 & year <= 2014]
month <- month[year >= 1950 & year <= 2014]
year <- year[year >= 1950 & year <= 2014]

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

# limit x
# SST <- SST[(x>=120 & x <= 250),,]

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

# identify columns in SST matrix corresponding to land
land <- is.na(colMeans(SST)) 

# For analysis, we only use the columns of the matrix with non-missing values:
X <- SST[,!land] 

# function to compute monthly means for a single time series
f <- function(x) tapply(x, month, mean, na.rm = TRUE)  

mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)

mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

anom <- X - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}


# get a vector of weights (square root of the cosine of latitude)
weight <- sqrt(cos(lat*pi/180))


# EOF
# weighting the columns
EOF.all <- svd.triplet(cov(anom.detr), col.w=weight) # this wasn't detrended in the last version of the script,
# might be a source of the problem?

# get % variance explained by era
var.all <- 100*round(prop.table(EOF.all$vs),3)

# get loadings for EOF1 and scale (reversing sign to match PDO)
eig.1.all <- -scale(EOF.all$U[,1])[,1]

# and get PC1 time series (reversing sign to match PDO)
pc1.all = -scale(anom.detr %*% EOF.all$U[,1])

# plot.dat <- cbind(eig.1.all, lon, lat) %>% ## COULD THIS CAUSE THE ERROR??
# as.data.frame(.) %>%
# mutate(eig.1.all * -1)

lim <- range(eig.1.all)

# plot loadings

png("./Figures/FCM_SST_EOF1_example.png", res = 300, width = 6, height = 5, units = 'in')

z <- rep(NA, ncol(SST))
z[!land] <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'South Korea', 'North Korea', 'Japan') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1950-2014 (", var.all[1], "%)", sep=""), cex=0.8)

dev.off()

# plot PC1 time series - this looks legit as it captures the 1976/77 PDO switch

# get decimal year and combine with PC1 in a DF to plot
plot.pc1 <- data.frame(decimal_year = year + (month - 0.5)/12,
                       PC1 = pc1.all,
                       "Rolling_13-month_mean" = rollmean(pc1.all, 13, align = "center", fill = NA)) %>%
  pivot_longer(cols = -decimal_year)

ggplot(plot.pc1, aes(decimal_year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = c("dark grey", "red")) 

ggsave("./Figures/FCM_SST_PC1_time_series_example.png", width = 6, height = 4, units = 'in')

## same thing with FCM SLP fields - examine EOF1/PC1 to make sure they look legit----------

# Extract time info for processing below (same across files)
files <- list.files(fcm.slp.dir, full.names = TRUE)

# Calculate EOF on SLP data ----
# first, load data
nc.slp <- nc_open(files[1]) # Isolate one CESM member

# this is clunky - get the correct time variable for CESM using ncdf4.helper
nc.slp <- nc_open(files[1])
d <- nc.get.time.series(nc.slp)

# names(dim(SLP)) <- c("lon", "lat", "time")

# now trim SLP to 1950:2014 
# Split each string at "-" and put results into a matrix
date_parts <- str_split_fixed(d, "-", n = 3)

# Extract year and month as separate vectors
year <- as.numeric(date_parts[, 1])
month <- as.numeric(date_parts[, 2])

# reload the NetCDF file
nc.slp <- tidync(files[1])

# Apply spatial filter on lat and lon using hyper_filter

# and also filter for 1950-2014

wanted <- which(year >= 1950 & year <= 2014)

# limit d, month, year to same years
d <- d[year >= 1950 & year <= 2014]
month <- month[year >= 1950 & year <= 2014]
year <- year[year >= 1950 & year <= 2014]

SLP <- nc.slp %>%
  hyper_filter(lat = between(lat, 20, 70),
               lon = between(lon, 120, 250),
               time = index %in% wanted) %>% 
  hyper_array(select_var = "PSL")

# Filter for the selected lat/lon after hyper_filter and pull their values as x,y

lon_tbl <- attr(SLP, "transforms")$lon

x <- lon_tbl %>% 
  filter(selected) %>% 
  pull(lon)


lat_tbl <- attr(SLP, "transforms")$lat

y <- lat_tbl %>% 
  filter(selected) %>% 
  pull(lat)

x; y

SLP <- SLP$PSL

dim(SLP) 
names(dim(SLP)) <- c("lon", "lat", "time")
# 105 lon x 53 lat x 780 months

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))


dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


# and plot 
SLP.mean <- colMeans(SLP)
z <- t(matrix(SLP.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(120,250))
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# that's totally right!

f <- function(x) tapply(x, month, mean, na.rm = TRUE)  # function to compute monthly means for a single time series

mu <- apply(SLP, 2, f)	# compute monthly means for each time series (cell)

mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

anom <- SLP - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}


# get a vector of weights (square root of the cosine of latitude)
weight <- sqrt(cos(lat*pi/180))

# EOF
# weighting the columns
EOF.all <- svd.triplet(cov(anom.detr), col.w=weight) # this wasn't detrended in the last version of the script,
# might be a source of the problem?

# get % variance explained by era
var.all <- 100*round(prop.table(EOF.all$vs),3)

# get loadings for EOF1 and scale 
eig.1.all <- scale(EOF.all$U[,1])[,1]

# and get PC1 time series 
pc1.all = scale(anom.detr %*% EOF.all$U[,1])

# plot.dat <- cbind(eig.1.all, lon, lat) %>% ## COULD THIS CAUSE THE ERROR??
# as.data.frame(.) %>%
# mutate(eig.1.all * -1)

lim <- range(eig.1.all)

# plot loadings - these are totally legit as they capture the AL

png("./Figures/FCM_SLP_EOF1_example.png", res = 300, width = 6, height = 5, units = 'in')

z <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'South Korea', 'North Korea', 'Japan') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1950-2014 (", var.all[1], "%)", sep=""), cex=0.8)

dev.off()

# plot PC1 time series - white noise, which is what we would expect

# get decimal year and combine with PC1 in a DF to plot
plot.pc1 <- data.frame(decimal_year = year + (month - 0.5)/12,
                       PC1 = pc1.all,
                       "Rolling_13-month_mean" = rollmean(pc1.all, 13, align = "center", fill = NA)) %>%
  pivot_longer(cols = -decimal_year)

ggplot(plot.pc1, aes(decimal_year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = c("dark grey", "red")) 

ggsave("./Figures/FCM_SLP_PC1_time_series_example.png", width = 6, height = 4, units = 'in')
