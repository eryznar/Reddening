# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# several packages that are helpful for this process have been removed from CRAN, 
# apparently because dependency PCICt was itself removed

# these downloads require package "remotes"

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

d <- nc.sst

# Load the NetCDF file
nc.sst <- tidync(files[1])

# Apply spatial filter on lat and lon using hyper_filter
SST <- nc.sst %>%
  hyper_filter(lat = between(lat, 20, 70),
               lon = between(lon, 120, 250)) %>% 
  hyper_array(select_var = "SST")

SST <- SST$SST

dim(SST) 
names(SST) <- c("lon", "lat", "time")
# 105 lon x 53 lat x 1980 months


# get x and y for plotting 
x <- dimnames(SST)[1]
y <- dimnames(SST)[2]

# wrangle for plotting
x <- as.numeric(x[[1]])
y <- as.numeric(y[[1]])

# this is clunky - get the correct time variable for CESM using ncdf4.helper
nc.sst <- nc_open(files[1])
d <- nc.get.time.series(nc.sst)

names(dim(SST)) <- c("lon", "lat", "time")

# now trim SST to 1950:2014 (when we expect EOF to be the PDO!)
# Split each string at "-" and put results into a matrix
date_parts <- str_split_fixed(d, "-", n = 3)

# Extract year and month as separate vectors
year <- as.numeric(date_parts[, 1])
month <- as.numeric(date_parts[, 2])

# limit SST to the desired years
SST <- SST[, , year %in% 1950:2014] 

# now get monthly means!
monthly_means <- Season(
  data = SST,
  time_dim = "time",   # specify your time dimension name
  monini = 1,           # first month of your data time series (e.g. January)
  moninf = 1,           # start month for averaging (January)
  monsup = 1,          # end month for averaging (December)
  method = mean,        # function to aggregate data
  na.rm = TRUE
)


dimnames(monthly_means)

dim(monthly_means)

str(monthly_means)

names(dim(monthly_means))

View(SST[,1,1])
# extract study area
# 20-70 deg. N, 120-250 deg. E
x <- ncvar_get(nc.sst, "lon")
y <- ncvar_get(nc.sst, "lat")

x; y # global, will need to be limited to study area

SST <- ncvar_get(nc.sst, "SST", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# limit x
SST <- SST[(x>=120 & x <= 250),,]

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
# SST.clean <- SST[, colSums(!is.na(SST)) > 0] # filtering columns with just NAs; COULD THIS CAUSE THE ERROR?

# identify columns in SST matrix corresponding to land
land <- is.na(colMeans(SST)) 

# For analysis, we only use the columns of the matrix with non-missing values:
X <- SST[,!land] 


f <- function(x) tapply(x, m, mean, na.rm = TRUE)  # function to compute monthly means for a single time series
# mu <- apply(SST.clean, 2, f)	# compute monthly means for each time series (cell)

mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)

mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 


# anom <- SST.clean - mu   # compute matrix of anomalies
anom <- X - mu   # compute matrix of anomalies

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

lim <- range(eig.1.all)

# Full time series!
z <- rep(NA, ncol(SST))
z[!land] <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Mexico') ,fill=T, add=T, lwd=1, col="darkgoldenrod3")
mtext(paste("EOF1 1949-2021 (", var.all[1], "%)", sep=""), cex=0.8)
text(230, 63, labels = paste("r (PDO) = ", round(cor(pcs$pc1.all, pcs$pdo), 2), sep = ""), cex = cex)


ggplot(plot.dat, aes(lon, lat, fill = eig.1.all))+
  geom_tile()+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgoldenrod", color = "black")+
  coord_cartesian(ylim = c(20, 67.4), xlim = c(125, 260), expand = FALSE)+
  xlab("Latitude")+
  ggtitle("FCM SSTa EOF1")+
  ylab("Longitude")+
  scale_fill_gradientn(colors = new.col, limits = c(-(max(na.omit(abs(plot.dat$eig.1.all)))), max(na.omit(abs(plot.dat$eig.1.all)))))




