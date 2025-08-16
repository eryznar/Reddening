#ERRST SST and NCEP/NCAR SLP processing

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

new.col <- oceColorsPalette(64)
### Process SST -----------------------------------------------------------------------------------------------------
# Load data
# nc.sst <- nc_open(paste0(dir, "Data/nceiErsstv5_ee08_74ee_6f8f.nc"))

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
  # ebs.x <- c(183, 183, 203, 203, 191)
  # #ebs.x <- ifelse(ebs.x > 180, ebs.x-360, ebs.x)
  # ebs.y <- c(53, 65, 65, 57.5, 53)
  
  # EBS polygon
  ebs.x <- c(187, 187, 203, 203, 191) 
  ebs.y <- c(53, 61, 61, 57.5, 53)

  
  xp <- cbind(ebs.x, ebs.y)
  loc=cbind(lon, lat)
  check <- in.poly(loc, xp=xp)
  
  SST.ebs[,!check] <- NA
  
  # and plot to check
  SST.mean <- colMeans(SST.ebs)
  z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(125,255))
  contour(x, y, z, add=T, col="white")  
  map('world2Hires',fill=F,add=T, lwd=2)
  
  
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
  
  SST.ebs <- SST.m[, 1:(ncol(SST.m)-2)]
  
  # and detrend
  SST.anom.ebs.detr <- SST.anom.ebs
  for(i in 1:ncol(SST.anom.ebs)) {
    xx = seq(1,nrow(SST.ebs))
    SST.anom.ebs.detr[,i] = SST.anom.ebs[,i] - predict(lm(SST.anom.ebs[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
  }
  
  
  # get average anomaly across the area
  SST.anom.ebs.detr <- rowMeans(SST.anom.ebs.detr)
  SST.ebs <- rowMeans(SST.ebs)
  
  # fit to winter means
  win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
  SST.win.anom.ebs.detr <- SST.anom.ebs.detr[m %in% c(11,12,1:3)]
  win.yr <- win.yr[m %in% c(11,12,1:3)]

  # also get winter absolute temps as a check
  SST.win.ebs <- SST.ebs[m %in% c(11,12,1:3)]
  
  plot.check <- data.frame(year = unique(win.yr),
                           ebs.win.mean = tapply(SST.win.ebs, win.yr, mean))  
  
  ggplot(plot.check, aes(year, ebs.win.mean)) +
    geom_line() # looks good
 
  
  
  data.frame(Date = names(SST.anom.ebs.detr), sst.anom = SST.anom.ebs.detr) %>%
    mutate(Year = as.numeric(as.character(substr(Date, 1, 4))),
           Month = as.numeric(as.character(substr(Date, 6, 7))),
           Win.year = case_when((Month %in% c(10:12)) ~ (Year+1),
                                TRUE ~ Year)) -> SST.anom.ebs.detr2
  
  rownames(SST.anom.ebs.detr2) <- NULL
  
  SST.anom.ebs.detr2%>%
    dplyr::select(!c(Date, Win.year)) %>%
    rename(month.anom = sst.anom) -> month.ebs.detr
  
  write.csv(month.ebs.detr, paste0("Output/ebs.monthlySSTanomalies_detrended.csv"))
  
  SST.anom.ebs.detr2 %>%
    group_by(Year) %>%
    reframe(mean.anom = mean(sst.anom)) -> SST.anom.ebs.detr.regyr
  
  # SST.anom.ebs.detr2 %>%
  #   group_by(Win.year) %>%
  #   reframe(mean.anom = mean(sst.anom)) %>%
  #   rename(Year = Win.year)-> SST.anom.ebs.detr.winyr
  
  
  SST.anom.ebs.detr2 %>%
    filter(Month %in% c(11:12, 1:3)) %>%
    group_by(Win.year) %>%
    reframe(mean.anom = mean(sst.anom)) %>%
    rename(Year = Win.year)-> SST.anom.ebs.detr.winter
  
  # write csv
  write.csv(SST.anom.ebs.detr.regyr, "./Output/SST.anom.ebs.detr.csv") # regular years
  # write.csv(SST.anom.ebs.detr.winyr, paste0(dir, "Output/SST.anom.ebs.winyr.detr.csv")) # Oct-Sept, year of January
  write.csv(SST.anom.ebs.detr.winter, "./Output/SST.winter.anom.ebs.detr.csv") # winter months, year of January
  
  ggplot(SST.anom.ebs.detr.winter, aes(Year, mean.anom)) + 
    geom_line()
  
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
  

  # and plot to check
  SST.mean <- colMeans(SST.goa)
  z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
  image(x,y,z, col=new.col, ylim=c(20,68), xlim=c(125,255))
  contour(x, y, z, add=T, col="white")  
  map('world2Hires',fill=F,add=T, lwd=2)
  
  
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
  
  SST.goa <- SST.m[, 1:(ncol(SST.m)-2)]
  
  # and detrend
  SST.anom.goa.detr <- SST.anom.goa
  for(i in 1:ncol(SST.anom.goa)) {
    xx = seq(1,nrow(SST.goa))
    SST.anom.goa.detr[,i] = SST.anom.goa[,i] - predict(lm(SST.anom.goa[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
  }
  
  
  
  # get average anomaly across the area
  SST.anom.goa.detr <- rowMeans(SST.anom.goa.detr)
  SST.goa <- rowMeans(SST.goa)
  
  # fit to winter means
  win.yr <- ifelse(m %in% c(11,12), yr+1, yr)
  SST.win.anom.goa.detr <- SST.anom.goa.detr[m %in% c(11,12,1:3)]
  win.yr <- win.yr[m %in% c(11,12,1:3)]

  # also get winter absolute temps as a check
  SST.win.goa <- SST.goa[m %in% c(11,12,1:3)]
  
  plot.check <- data.frame(year = unique(win.yr),
                           goa.win.mean = tapply(SST.win.goa, win.yr, mean))  
  
  ggplot(plot.check, aes(year, goa.win.mean)) +
    geom_line() # looks good
  
  
  data.frame(Date = names(SST.anom.goa.detr), sst.anom = SST.anom.goa.detr) %>%
    mutate(Year = as.numeric(as.character(substr(Date, 1, 4))),
           Month = as.numeric(as.character(substr(Date, 6, 7))),
           Win.year = case_when((Month %in% c(10:12)) ~ (Year+1),
                                TRUE ~ Year)) -> SST.anom.goa.detr2
  
  rownames(SST.anom.goa.detr2) <- NULL
  
  SST.anom.goa.detr2%>%
    dplyr::select(!c(Date, Win.year)) %>%
    rename(month.anom = sst.anom) -> month.goa.detr
  
  write.csv(month.goa.detr, paste0("Output/goa.monthlySSTanomalies_detrended.csv"))
  
  SST.anom.goa.detr2 %>%
    group_by(Year) %>%
    reframe(mean.anom = mean(sst.anom)) -> SST.anom.goa.detr.regyr
  
  # SST.anom.goa.detr2 %>%
  #   group_by(Win.year) %>%
  #   reframe(mean.anom = mean(sst.anom)) %>%
  #   rename(Year = Win.year)-> SST.anom.goa.detr.winyr
  
  
  SST.anom.goa.detr2 %>%
    filter(Month %in% c(11:12, 1:3)) %>%
    group_by(Win.year) %>%
    reframe(mean.anom = mean(sst.anom)) %>%
    rename(Year = Win.year)-> SST.anom.goa.detr.winter
  
  # write csv
  write.csv(SST.anom.goa.detr.regyr, "./Output/SST.anom.goa.detr.csv") # regular years
  # write.csv(SST.anom.goa.detr.winyr, "./Output/SST.anom.goa.winyr.detr.csv") # Oct-Sept, year of January
  write.csv(SST.anom.goa.detr.winter, "./Output/SST.winter.anom.goa.detr.csv") # winter months, year of January
  

### Process SLP -----------------------------------------------------------------------------------------------------
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

# save monthly anomalies
data.frame(Date = names(SLP.anom), slp = SLP.anom) %>%
  mutate(Year = as.numeric(as.character(substr(Date, 1, 4))),
         Month = as.numeric(as.character(substr(Date, 6, 7))),
         Win.year = case_when((Month %in% c(10:12)) ~ (Year+1),
                              TRUE ~ Year)) -> SLP.anom2

rownames(SLP.anom2) <-NULL

SLP.anom2%>%
  dplyr::select(!c(Date, Win.year)) %>%
  rename(month.anom = slp) -> month.slp

write.csv(month.slp, paste0(dir, "Output/monthlySLPanomalies.csv"))


# Calculate winter means
SLP.anom2 %>%
  filter(Month %in% c(11, 12, 1:3)) %>%
  group_by(Win.year) %>%
  reframe(SLP.win.anom = mean(slp)) -> SLP.dat


write.csv(SLP.dat, paste0(dir, "Output/monthlywinterSLPanomalies.csv"))



### Calculate EOF on SLP data ----
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

slp.anom <- SLP.m[,1:35] - mu 


# get a vector of weights (square root of the cosine of latitude)

# identify columns containing NA
temp1 <- str_split(colnames(SLP.m), "E", simplify = T)[,1]
X.lats <- as.numeric(str_split(temp1, "N", simplify = T)[,2])

weight <- sqrt(cos(X.lats*pi/180)) 


pca <- FactoMineR::svd.triplet(cov(slp.anom), col.w=na.omit(weight)) #weighting the columns

pc1_slp <- as.matrix(slp.anom) %*% pca$U[,1]

# and scale!
pc1_slp <- as.vector(scale(pc1_slp))

write.csv(pc1_slp, paste0(dir, "Output/PC1slp.monthlyanomalies.csv"))

