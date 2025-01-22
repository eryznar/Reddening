### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

source("./Scripts/load.libs.functions.R")
source("Y:/KOD_Survey/EBS Shelf/Spatial crab/load.spatialdata.R")

### READ IN WINTER SSTa ----------------------------------
ebs.SSTa <- read.csv(paste0(dir, "Output/SST.winter.anom.ebs.csv")) %>%
  dplyr::select(!X)
goa.SSTa <- read.csv(paste0(dir, "Output/SST.winter.anom.goa.csv")) %>%
  dplyr::select(!X)


### PROCESS SLPa for NORTH PACIFIC REGION -----------------------------------------
# first, load data
nc.slp <- nc_open(paste0(dir, "Data/hawaii_soest_f19d_3925_d70b_1f05_6ec6_fbce5-90N.nc"))

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
poly.x <- c(130, 130, 250, 250, 130) 
poly.y <- c(5, 70, 70, 5, 5)

xp <- cbind(poly.x, poly.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

SLP[,!check] <- NA

# plot to check
z <- colMeans(SLP*100)
z <- t(matrix(z, length(y)))
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", ylim=c(5,70), xlim=c(130,250))
contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires',fill=F, xlim=c(130,250), ylim=c(5,70),add=T, lwd=1)
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

slp.anom <- cbind(SLP.m[,1:(ncol(SLP.m)-2)]  - mu, data.frame(year = SLP.m$year, month = SLP.m$month))

# Pivot longer into a format I recognize!
SLPa.long <- pivot_longer(slp.anom, cols = c(-month, -year), names_to = "crds", values_to = "anom") %>%
                mutate(lat = str_split(crds, "E", simplify = T)[,1],
                       lat = as.numeric(str_split(lat, "N", simplify = T)[,2]), # pull out lat
                       lon = as.numeric(str_split(crds, "E", simplify = T)[,2]), #pull out lon
                       #lon = ifelse(lon > 180, lon-360, lon),# convert lon
                       month = as.numeric(month),
                       year = as.numeric(year)) %>% 
                filter(month %in% c(11:12, 1:3)) %>% # filter to winter months
                group_by(crds, lat, lon, year) %>%
                reframe(slp.win.anom = mean(anom)) 
            


# Filter by strong and weak red periods
SLPa.red <- SLPa.long %>%
  filter(year %in% c(2005:2024)) 

SLPa.white <- SLPa.long %>%
  filter(year %in% c(1975:1995)) 

# Get map layers
mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))



# Plot to check (does this look correct??)
SLPa.red2 <- SLPa.red %>%
              group_by(lat, lon) %>%
              reframe(slp.win.anom = mean(slp.win.anom))

SLPa.white2 <- SLPa.white %>%
  group_by(lat, lon) %>%
  reframe(slp.win.anom = mean(slp.win.anom))


ggplot() +
  geom_tile(SLPa.red, mapping = aes(lon, lat, fill = slp.win.anom))+
  stat_contour(SLPa.red, mapping = aes(lon, lat, z = slp.win.anom ), geom = "contour", position = "identity")+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(15, 70), xlim = c(160, 250), expand = FALSE)+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
  ggtitle("Winter SLPa during red noise (2005-2024)")

ggplot() +
  geom_tile(SLPa.white2, mapping = aes(lon, lat, fill = slp.win.anom))+
  stat_contour(SLPa.white2, mapping = aes(lon, lat, z = slp.win.anom ), geom = "contour", position = "identity")+
  geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "lightgrey", color = "darkgrey")+
  coord_sf(ylim = c(15, 70), xlim = c(160, 250), expand = FALSE)+
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
  ggtitle("Winter SLPa during white noise (1975-1995)")

# Regressions
  #GOA
  goa.red <- data.frame()
  goa.white <- data.frame()
  
  sst <- goa.SSTa
  
  identical(unique(SLPa.red$crds), unique(SLPa.white$crds))
  
  cc <- unique(SLPa.red$crds)
  
  for(ii in 1:length(cc)){
    SLP.r <- SLPa.red %>% filter(crds== cc[ii])
    SLP.w <- SLPa.white %>% filter(crds== cc[ii])
    
    sst.r <- sst %>% filter(Year %in% SLP.r$year)
    sst.w <- sst %>% filter(Year %in% SLP.w$year)
    
    # Red noise
    mod.r <- lm(SLP.r$slp.win.anom ~ sst.r$mean.sst)
    est <- summary(mod.r)$coef[2,1]
    r <- sqrt(summary(mod.r)$r.squared)
    
    df.r <- data.frame(crds = SLP.r$crds, 
                       lat = SLP.r$lat, 
                       lon = SLP.r$lon,
                       est = est,
                       r= r)
    
    goa.red <- rbind(goa.red, df.r)
    
    # White noise
    mod.w <- lm(SLP.w$slp.win.anom ~ sst.w$mean.sst)
    est <- summary(mod.w)$coef[2,1]
    r <- sqrt(summary(mod.w)$r.squared)
    
    df.w <- data.frame(crds = SLP.w$crds, 
                       lat = SLP.w$lat, 
                       lon = SLP.w$lon,
                       est = est,
                       r = r)
    
    goa.white <- rbind(goa.white, df.w)
    
  }
  
  # Plot regressions
  ggplot() +
    geom_tile(goa.red, mapping = aes(lon, lat, fill = est))+
    stat_contour(goa.red, mapping = aes(lon, lat, z = est ), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("GOA SLP vs. SST (red noise: 2005-2024)")+
    theme_bw()-> r1
  
  ggplot() +
    geom_tile(goa.white, mapping = aes(lon, lat, fill = est))+
    stat_contour(goa.white, mapping = aes(lon, lat, z = est ), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("GOA SLP vs. SST (white noise: 1975-1995)")+
    theme_bw() -> r2
  
  ggplot() +
    geom_tile(goa.red, mapping = aes(lon, lat, fill = r))+
    stat_contour(goa.red, mapping = aes(lon, lat, z = r), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("GOA SLP vs. SST (red noise: 2005-2024)") +
    theme_bw() -> r3
  
  ggplot() +
    geom_tile(goa.white, mapping = aes(lon, lat, fill = r))+
    stat_contour(goa.white, mapping = aes(lon, lat, z = r), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("GOA SLP vs. SST (white noise: 1975-1995)")+
    theme_bw() -> r4
  
 
  r1 + r3 + r2 + r4 +  plot_layout(ncol = 2) #warnings are ok because the contour is the same
  
  ggsave("./Figures/goaSLPv.SSTregression.png", width = 8, height = 5, units = "in")
  
  # EBS
  ebs.red <- data.frame()
  ebs.white <- data.frame()
  
  sst <- ebs.SSTa
  
  identical(unique(SLPa.red$crds), unique(SLPa.white$crds))
  
  cc <- unique(SLPa.red$crds)
  
  for(ii in 1:length(cc)){
    SLP.r <- SLPa.red %>% filter(crds== cc[ii])
    SLP.w <- SLPa.white %>% filter(crds== cc[ii])
    
    sst.r <- sst %>% filter(Year %in% SLP.r$year)
    sst.w <- sst %>% filter(Year %in% SLP.w$year)
    
    # Red noise
    mod.r <- lm(SLP.r$slp.win.anom ~ sst.r$mean.sst)
    est <- summary(mod.r)$coef[2,1]
    r <- sqrt(summary(mod.r)$r.squared)
    
    df.r <- data.frame(crds = SLP.r$crds, 
                       lat = SLP.r$lat, 
                       lon = SLP.r$lon,
                       est = est,
                       r = r)
    
    ebs.red <- rbind(ebs.red, df.r)
    
    # White noise
    mod.w <- lm(SLP.w$slp.win.anom ~ sst.w$mean.sst)
    est <- summary(mod.w)$coef[2,1]
    r <- sqrt(summary(mod.w)$r.squared)
    
    df.w <- data.frame(crds = SLP.w$crds, 
                       lat = SLP.w$lat, 
                       lon = SLP.w$lon,
                       est = est,
                      r = r)
    
    ebs.white <- rbind(ebs.white, df.w)
    
  }
  
  # Plot regressions
  ggplot() +
    geom_tile(ebs.red, mapping = aes(lon, lat, fill = est))+
    stat_contour(ebs.red, mapping = aes(lon, lat, z = est ), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("EBS SLP vs. SST (red noise: 2005-2024)")+
    theme_bw()-> r1
  
  ggplot() +
    geom_tile(ebs.white, mapping = aes(lon, lat, fill = est))+
    stat_contour(ebs.white, mapping = aes(lon, lat, z = est ), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("EBS SLP vs. SST (white noise: 1975-1995)")+
    theme_bw() -> r2
  
  ggplot() +
    geom_tile(ebs.red, mapping = aes(lon, lat, fill = r))+
    stat_contour(ebs.red, mapping = aes(lon, lat, z = r), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("EBS SLP vs. SST (red noise: 2005-2024)") +
    theme_bw() -> r3
  
  ggplot() +
    geom_tile(ebs.white, mapping = aes(lon, lat, fill = r))+
    stat_contour(ebs.white, mapping = aes(lon, lat, z = r), geom = "contour", position = "identity", color = "lightgrey")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = "darkgrey", color = "black")+
    coord_sf(ylim = c(27, 70), xlim = c(160, 250), expand = FALSE)+
    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    ggtitle("EBS SLP vs. SST (white noise: 1975-1995)")+
    theme_bw() -> r4
  
  
  r1 + r3 + r2 + r4 +  plot_layout(ncol = 2) #warnings are ok because the contour is the same
  
  ggsave("./Figures/ebsSLPv.SSTregression.png", width = 8, height = 5, units = "in")
  
