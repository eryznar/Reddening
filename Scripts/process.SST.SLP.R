#ERRST SST and NCEP/NCAR SLP processing

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

### SST -----------------------------------------------------------------------------------------------------
# Load data
nc.sst <- nc_open("./Data/nceiErsstv5_ee08_74ee_6f8f.nc")

# process sst data - first, extract dates
raw <- ncvar_get(nc.sst, "time")  # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract coordinates
x <- ncvar_get(nc.sst, "longitude")
y <- ncvar_get(nc.sst, "latitude")

# extract data
SST <- ncvar_get(nc.sst, "sst", verbose = F)

# Change data to a matrix
SST <- aperm(SST, 3:1)  

# Filter to area of interest
  # EBS: ----
  # Change to matrix
  SST.ebs <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  
  
  # Get lat/long vectors and add names to SST matrix
  lat <- rep(y, length(x))   
  lon <- rep(x, each = length(y))   
  dimnames(SST.ebs) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
  
  # Filter to region
  ebs.x <- c(183, 183, 203, 203, 191)
  #ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
  ebs.y <- c(53, 65, 65, 57.5, 53)
  
  xp <- cbind(ebs.x, ebs.y)
  loc=cbind(lon, lat)
  check <- in.poly(loc, xp=xp)
  
  SST.ebs[,!check] <- NA
  
  # plot to check
  z <- colMeans(SST.ebs)
  z <- t(matrix(z, length(y)))
  image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", ylim=c(50,66), xlim=c(180,208))
  contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
  map('world2Hires',fill=F, ylim=c(50,66), xlim=c(180,208),add=T, lwd=1)
  # looks good
  
  SST.ebs %>%
    as.data.frame(.) %>%
    mutate(date = rownames(.)) %>%
    pivot_longer(!date, names_to = "coords", values_to = "SST") %>%
    na.omit() %>%
    pivot_wider(., names_from = coords, values_from = SST) %>%
    as.data.frame(.) -> clean.SST.ebs
  
  rownames(clean.SST.ebs) = clean.SST.ebs$date
  
  clean.SST.ebs <- clean.SST.ebs[,-1]
  
  # now we need to get monthly means!
  m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates
  
  f <- function(x) tapply(x, m.y, mean)
  SST.m <- as.data.frame(apply(clean.SST.ebs, 2, f))
  
  
  vv <- matrix(unlist(strsplit(as.character(rownames(SST.m)), "-")),ncol=2, byrow = T)
  
  SST.m$year <- as.numeric(vv[,1])
  SST.m$month <- as.numeric(vv[,2])
  
  SST.m <- SST.m %>%
    arrange(month) %>%
    arrange(year)
  
  # remove seasonal signal
  m <- SST.m$month
  yr <- SST.m$year
  
  f <- function(x) tapply(x, m, mean)
  mu <- apply(SST.m, 2, f)	# Compute monthly means for each cell
  
  # process as for SST
  mu <- mu[rep(1:12, floor(length(d)/12)),] 
  xtra <- 12*((length(d)/12)-floor(length(d)/12))
  mu <- rbind(mu, mu[1:xtra,])
  
  SST.anom.ebs <- SST.m[,1:(ncol(SST.m)-2)] - mu   # Compute matrix of anomalies - dropping year and month!
  
  # get average anomaly across the area
  SST.anom.ebs <- rowMeans(SST.anom.ebs)
  
  # fit to winter means
  win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
  SST.win.anom.ebs <- SST.anom.ebs[m %in% c(11,12,1:3)]
  win.yr <- win.yr[m %in% c(11,12,1:3)]
  SST.win.anom.ebs <- tapply(SST.win.anom.ebs, win.yr, mean)
  
  plot(1854:2024, SST.win.anom.ebs, type = "l")
  
  # pivot longer
  data.frame(Year = names(SST.win.anom.ebs), mean.sst = SST.win.anom.ebs) -> SST.win.anom.ebs
  
  rownames(SST.win.anom.ebs) <- NULL
  
  # write csv
  write.csv(SST.win.anom.ebs, "./Output/SST.winter.anom.ebs.csv")
  
  
  # GOA: ----
  # Change to matrix
  SST.goa <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  
  
  # Get lat/long vectors and add names to SST matrix
  lat <- rep(y, length(x))   
  lon <- rep(x, each = length(y))   
  dimnames(SST.goa) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))
  
  # Filter to region
  goa.x <- c(201, 201, 205, 208, 225, 231, 201)
  #goa.x <- ifelse(goa.x > 180, goa.x-360, goa.x)
  goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)
  
  xp <- cbind(goa.x, goa.y)
  loc=cbind(lon, lat)
  check <- in.poly(loc, xp=xp)
  
  SST.goa[,!check] <- NA
  
  # plot to check
  z <- colMeans(SST.goa)
  z <- t(matrix(z, length(y)))
  image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", ylim=c(53,66), xlim=c(199,230))
  contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
  map('world2Hires',fill=F, ylim=c(53,66), xlim=c(199,230),add=T, lwd=1)
  # looks good
  
  SST.goa %>%
    as.data.frame(.) %>%
    mutate(date = rownames(.)) %>%
    pivot_longer(!date, names_to = "coords", values_to = "SST") %>%
    na.omit() %>%
    pivot_wider(., names_from = coords, values_from = SST) %>%
    as.data.frame(.) -> clean.SST.goa
  
  rownames(clean.SST.goa) = clean.SST.goa$date
  
  clean.SST.goa <- clean.SST.goa[,-1]
  
  # now we need to get monthly means!
  m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates
  
  f <- function(x) tapply(x, m.y, mean)
  SST.m <- as.data.frame(apply(clean.SST.goa, 2, f))
  
  
  vv <- matrix(unlist(strsplit(as.character(rownames(SST.m)), "-")),ncol=2, byrow = T)
  
  SST.m$year <- as.numeric(vv[,1])
  SST.m$month <- as.numeric(vv[,2])
  
  SST.m <- SST.m %>%
    arrange(month) %>%
    arrange(year)
  
  # remove seasonal signal
  m <- SST.m$month
  yr <- SST.m$year
  
  f <- function(x) tapply(x, m, mean)
  mu <- apply(SST.m, 2, f)	# Compute monthly means for each cell
  
  # process as for SST
  mu <- mu[rep(1:12, floor(length(d)/12)),] 
  xtra <- 12*((length(d)/12)-floor(length(d)/12))
  mu <- rbind(mu, mu[1:xtra,])
  
  SST.anom.goa <- SST.m[,1:(ncol(SST.m)-2)] - mu   # Compute matrix of anomalies - dropping year and month!
  
  # get average anomaly across the area
  SST.anom.goa <- rowMeans(SST.anom.goa)
  
  # fit to winter means
  win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
  SST.win.anom.goa <- SST.anom.goa[m %in% c(11,12,1:3)]
  win.yr <- win.yr[m %in% c(11,12,1:3)]
  SST.win.anom.goa <- tapply(SST.win.anom.goa, win.yr, mean)
  
  # Plot
  plot(1854:2024, SST.win.anom.goa, type = "l")
  
  # pivot longer
  data.frame(Year = names(SST.win.anom.goa), mean.sst = SST.win.anom.goa) -> SST.win.anom.goa
  
  rownames(SST.win.anom.goa) <- NULL
  
  # write csv
  write.csv(SST.win.anom.goa, "./Output/SST.winter.anom.goa.csv")
  

### SLP -----------------------------------------------------------------------------------------------------
# first, load data
nc.slp <- nc_open("Data/hawaii_soest_f19d_3925_d70b_1322_e90d_09e0NEW.nc")

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
poly.x <- c(191, 191, 208, 208, 191) 
poly.y <- c(44, 55, 55, 44, 44)

xp <- cbind(poly.x, poly.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

SLP[,!check] <- NA

# plot to check
z <- colMeans(SLP*100)
z <- t(matrix(z, length(y)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", ylim=c(35,66), xlim=c(170,220))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(20,66),add=T, lwd=1)
# looks good

SLP %>%
  as.data.frame(.) %>%
  mutate(date = rownames(.)) %>%
  #mutate(date = mdy(rownames(.))) %>%
  pivot_longer(!date, names_to = "coords", values_to = "SLP") %>%
  # mutate(year = year(date),
  #        month = month(date)) %>%
  # mutate(year = case_when((year > 2024) ~ year - 100,
  #                        TRUE ~ year)) %>%
  na.omit() %>%
  mutate(SLP = SLP * 100) %>%
  pivot_wider(., names_from = coords, values_from = SLP) %>%
  as.data.frame(.) -> clean.SLP

rownames(clean.SLP) = clean.SLP$date

clean.SLP <- clean.SLP[,-1]

# now we need to get monthly means!
m.y <- paste(years(d), as.numeric(months(d)), sep="-") # make a month-year factor from the dates

f <- function(x) tapply(x, m.y, mean)
SLP.m <- as.data.frame(apply(clean.SLP, 2, f))

vv <- matrix(unlist(strsplit(as.character(rownames(SLP.m)), "-")),ncol=2, byrow = T)

SLP.m$year <- as.numeric(vv[,1])
SLP.m$month <- as.numeric(vv[,2])

SLP.m <- SLP.m %>%
  arrange(month) %>%
  arrange(year)

# remove seasonal signal
m <- SLP.m$month
yr <- SLP.m$year

f <- function(x) tapply(x, m, mean)
mu <- apply(SLP.m, 2, f)	# Compute monthly means for each cell

# process as for SST
mu <- mu[rep(1:12, floor(length(d)/12)),] 
xtra <- 12*((length(d)/12)-floor(length(d)/12))
mu <- rbind(mu, mu[1:xtra,])

SLP.anom <- SLP.m[,1:35] - mu   # Compute matrix of anomalies - dropping year and month!

# get average anomaly across the area
SLP.anom <- rowMeans(SLP.anom)

# fit to winter means
win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
SLP.win.anom <- SLP.anom[m %in% c(11,12,1:3)]
win.yr <- win.yr[m %in% c(11,12,1:3)]
SLP.win.anom <- tapply(SLP.win.anom, win.yr, mean)

#SLP.win.anom <- SLP.win.anom[2:72]
plot(1948:2024, SLP.win.anom, type="l")

# and calculate standard deviation over 11-year rolling windows
SLP.win.sd <- rollapply(SLP.win.anom, 11, sd, fill=NA)
plot(1948:2024, SLP.win.sd, type="l")

# now fit a non-parametric regression

# first, make a data frame
plot.dat <- data.frame(year= 1948:2024, sd=SLP.win.sd) %>%
  na.omit(.)

# fit the model
mod <- gam(sd ~ s(year), data=plot.dat)
pred <- predict(mod, se=T, newdata = plot.dat)
plot.dat$mean <- pred$fit                    

b.plot <- ggplot(plot.dat, aes(year, sd)) +
  geom_line(size=0.2) +
  geom_line(aes(year, mean), color="salmon", size=0.4) + # geom_ribbon(aes(ymin=LCI, ymax=UCI), alpha=0.2) +
  theme(axis.title.x = element_blank(), plot.title = element_text(size=8), axis.text = element_text(size=7),
        axis.title.y = element_text(size=7)) +
  ylab("Standard deviation (pa)") +
  ggtitle("Aleutian Low variability") +
  geom_vline(xintercept = 1988.5, lty=2, size=0.3) +
  xlim(1950,2024)

# save data
data.frame(year = names(SLP.win.anom), SLP.win.anom = SLP.win.anom) -> SLP.dat

rownames(SLP.dat) = NULL

write.csv(SLP.dat, "./Output/SLP.winter.anom.csv")
