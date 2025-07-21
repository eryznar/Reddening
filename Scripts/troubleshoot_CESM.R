# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# Get map layers
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))

yrs <- 1850:2013 # these are years that are similar across SLP and SST models

# Color palette
library(oce)

new.col <- oceColorsPalette(64)

# Identify folders to download files from
fcm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/FCM/") #FCM SST
mdm.sst.dir <- paste0(dir, "Data/CESM2 ensemble/SST/MDM/") #MDM SST

fcm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/FCM/") #FCM SLP
mdm.slp.dir <- paste0(dir, "Data/CESM2 ensemble/SLP/MDM/") #MDM SLP

# Extract time info for processing below (same across files)
files <- list.files(fcm.sst.dir, full.names = TRUE)

### Calculate EOF on SST data ----
# first, load data
nc.sst <- nc_open(files[1]) # Isolate one CESM member

# process SST data - first, extract dates (different from legacy script)
time_units <- ncmeta::nc_atts(files[1], "time") %>% 
  filter(name == "units") %>%
  pull(value)

unit_parts <- str_split(time_units, " since ")[[1]]
time_unit <- unit_parts[1]
origin_date <- ymd_hms(unit_parts[2])

raw <- ncvar_get(nc.sst, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
days <- lubridate::days(raw) + origin_date
d <- date(days)


# extract study area
# 20-70 deg. N, 120-250 deg. E
x <- ncvar_get(nc.sst, "lon")
y <- ncvar_get(nc.sst, "lat")

x; y


SST <- ncvar_get(nc.sst, "SST", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


# Filter to region
ebs.x <- c(125, 125, 255, 255, 125)
#ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
ebs.y <- c(20, 68, 68, 20, 20)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

SST[,!check] <- NA

m <- months(d)  # Extracts months from the date vector
yr <- years(d)
#m <- match(m, month.name)

# # reset lat/lon
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y))   

# and plot 
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(125,255))
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# remove seasonal means
SST.clean <- SST[, colSums(!is.na(SST)) > 0] # filtering columns with just NAs; COULD THIS CAUSE THE ERROR?

f <- function(x) tapply(x, m, mean, na.rm = TRUE)  # function to compute monthly means for a single time series
mu <- apply(SST.clean, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 


anom <- SST.clean - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}


# get a vector of weights (square root of the cosine of latitude)
temp1 <- str_split(colnames(anom), "E", simplify = T)[,1]
lat <- as.numeric(str_split(temp1, "N", simplify = T)[,2])
lon <- as.numeric(str_split(colnames(anom), "E", simplify = T)[,2])
weight <- sqrt(cos(lat*pi/180))


# EOF by era 
# weighting the columns
EOF.all <- svd.triplet(cov(anom), col.w=weight)

# get loadings for EOF1 and scale
eig.1.all <- scale(EOF.all$U[,1])[,1]

plot.dat <- cbind(eig.1.all, lon, lat) %>% ## COULD THIS CAUSE THE ERROR??
  as.data.frame(.) %>%
  mutate(eig.1.all * -1)



ggplot(plot.dat, aes(lon, lat, fill = eig.1.all))+
  geom_tile()+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
  xlab("Latitude")+
  ggtitle("FCM SSTa EOF1")+
  ylab("Longitude")+
  scale_fill_gradientn(colors = new.col, limits = c(-(max(na.omit(abs(plot.dat$eig.1.all)))), max(na.omit(abs(plot.dat$eig.1.all)))))




