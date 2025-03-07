library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(FactoMineR)
library(nlme)
library(MuMIn)
library(oce)
library(ggpubr)
source("./Scripts/load.libs.functions.R")

# load NCEP NCAR data

# source:
# https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_f19d_3925_d70b.nc?slp[(1948-01-01):1:(2024-12-01T00:00:00Z)][(20):1:(67.5)][(130):1:(250)]

# set palettes
new.col <- oceColorsPalette(64)
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set theme
theme_set(theme_bw())

# load and process SLP data
nc <- nc_open("./data/hawaii_soest_f19d_3925_d70b_53f6_4d1c_e436.nc")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 20-70 deg. N, 120-250 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

x; y


SLP <- ncvar_get(nc, "slp", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SLP <- aperm(SLP, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SLP <- matrix(SLP, nrow=dim(SLP)[1], ncol=prod(dim(SLP)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SLP) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


m <- months(d)  # Extracts months from the date vector
yr <- years(d)


# # reset lat/lon
# lat <- rep(y, length(x))   
# lon <- rep(x, each = length(y))   


# and plot 
SLP.mean <- colMeans(SLP)
z <- t(matrix(SLP.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col)
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# change to mb
# SLP <- SLP/100

# Plot annual winter average SLP fields 1948-2024
  SLP.long <- SLP %>%
    as.data.frame(.) %>%
    mutate(date = rownames(.),
           month = m,
           year = as.numeric(as.character(yr))) %>%
    pivot_longer(., cols=(!c("date", "month", "year")), names_to = "crds", values_to = "slp") %>%
    mutate(lat = str_split(crds, "E", simplify = T)[,1],
           lat = as.numeric(str_split(lat, "N", simplify = T)[,2]), # pull out lat
           lon = as.numeric(str_split(crds, "E", simplify = T)[,2])) %>% 
    filter(month %in% c("Nov", "Dec", "Jan", "Feb", "Mar")) %>%  #filter to winter months
    mutate(win.year = case_when((month %in% c("Nov", "Dec")) ~ year+1, 
                                TRUE ~ year)) %>% #specify winter year
    group_by(crds, lat, lon, win.year) %>%
    reframe(mean.slp.win = mean(slp)) 
  
  # Get map
  mapWorld <- map_data('world', wrap=c(-25,335), ylim=c(-55,75))
  
  
  # Plot 1989-2024
  SLP.long %>%
    filter(win.year %in% c(1989:2024)) -> SLP.long2
  
  ggplot() +
    geom_tile(SLP.long2, mapping = aes(lon, lat, fill = as.numeric(mean.slp.win)))+
    #scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    scale_fill_gradientn(colours = c(scales::muted("blue"),"white",scales::muted("red")), 
                         values = rescale(c(min(SLP.long$mean.slp.win),median(SLP.long$mean.slp.win),max(SLP.long$mean.slp.win))),
                         guide = "colorbar",
                         name = "Mean slp")+
    facet_wrap(~win.year)+
    geom_contour(SLP.long2, mapping = aes(lon, lat, z = mean.slp.win), position = "identity", color = "white")+
    # geomtextpath::geom_textcontour(SLP.long, mapping = aes(lon, lat, z = mean.slp.win), 
    #                                padding = unit(0.05, "in"),
    #                                color = "white")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = NA, color = "black")+
    coord_sf(ylim = c(20, 69), xlim = c(130, 250), expand = FALSE)+
    scale_y_continuous(breaks = seq(min(SLP.long$lat), max(SLP.long$lat), by = 15))+
    scale_x_continuous(breaks = seq(min(SLP.long$lat), max(SLP.long$lon), by = 40))+
    theme_bw()+
    ylab("Latitude") +
    xlab("Longitude")+
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 8)) -> mean.slp
  
  ggsave(plot = mean.slp, "./Figures/meanwinterSLPfields1989-2024.png", width = 8.5, height = 6, units = "in")

  # Plot 1948-2024
  SLP.long %>%
    filter(win.year %in% c(1948:2024)) -> SLP.long2
  
  ggplot() +
    geom_tile(SLP.long2, mapping = aes(lon, lat, fill = as.numeric(mean.slp.win)))+
    #scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = "white")+
    scale_fill_gradientn(colours = c(scales::muted("blue"),"white",scales::muted("red")), 
                         values = rescale(c(min(SLP.long$mean.slp.win),median(SLP.long$mean.slp.win),max(SLP.long$mean.slp.win))),
                         guide = "colorbar",
                         name = "Mean slp")+
    facet_wrap(~win.year)+
    geom_contour(SLP.long2, mapping = aes(lon, lat, z = mean.slp.win), position = "identity", color = "white")+
    # geomtextpath::geom_textcontour(SLP.long, mapping = aes(lon, lat, z = mean.slp.win), 
    #                                padding = unit(0.05, "in"),
    #                                color = "white")+
    geom_polygon(data = mapWorld, aes(x=long, y = lat, group = group), fill = NA, color = "black")+
    coord_sf(ylim = c(20, 69), xlim = c(130, 250), expand = FALSE)+
    scale_y_continuous(breaks = seq(min(SLP.long$lat), max(SLP.long$lat), by = 15))+
    scale_x_continuous(breaks = seq(min(SLP.long$lat), max(SLP.long$lon), by = 40))+
    theme_bw()+
    ylab("Latitude") +
    xlab("Longitude")+
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          strip.text = element_text(size = 8)) -> mean.slp
  
  ggsave(plot = mean.slp, "./Figures/meanwinterSLPfields1948-2024.png", width = 8.5, height = 7, units = "in")

# looks good!

# remove seasonal means
f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(SLP, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 


anom <- SLP - mu   # compute matrix of anomalies

# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}

# get a vector of weights (square root of the cosine of latitude)
weight <- sqrt(cos(lat*pi/180))

# restrict to NDJFM
SLP.NDJFM <- anom.detr[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec"),]

SLP.raw <- SLP[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec"),] # these just need to be assigned to winter year
# corresponding to January, then mean values can be plotted by year for 1989:2024 

# EOF by era 
# weighting the columns
EOF.all <- svd.triplet(cov(SLP.NDJFM), col.w=weight)



# get loadings for EOF1-2 by era
eig.1.all <- EOF.all$U[,1]
eig.2.all <- EOF.all$U[,2]



# now fit to white-noise period (1989:2004) and red-noise period (2005-2024)
yr.NDJFM <- yr[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec")]
  
EOF.1 <- svd.triplet(cov(SLP.NDJFM[yr.NDJFM %in% 1989:2004,]), col.w=weight)
EOF.2 <- svd.triplet(cov(SLP.NDJFM[yr.NDJFM %in% 2005:2024,]), col.w=weight)

# get loadings for EOF1-2 by era
eig.1.1 <- EOF.1$U[,1] 
eig.2.1 <- EOF.1$U[,2]

eig.1.2 <- EOF.2$U[,1] 
eig.2.2 <- EOF.2$U[,2]

# get % variance explained by era
var.all <- 100*round(prop.table(EOF.all$vs),3)
var.1 <- 100*round(prop.table(EOF.1$vs),3)
var.2 <- 100*round(prop.table(EOF.2$vs),3)

# set the limit for plotting 
lim.1 <- range(eig.1.all, eig.1.1, eig.1.2)
lim.2 <- range(eig.2.all, eig.2.1, eig.2.2)

lim.1
lim.2

#################
# Full time series!

png("./Figures/slp_EOF.png", 6, 3, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(1,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))


z  <- eig.1.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1948-2024 (", var.all[1], "%)", sep=""), cex=0.8)
polygon(x = c(186.25, 208.75, 208.75, 186.25), y = c(56.25, 56.25, 46, 46), border = "white", lwd = 2)

z  <- eig.2.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim.2[1], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1948-2024 (", var.all[2], "%)", sep=""), cex=0.8)

dev.off()

# Era 1!
z  <- eig.1.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 1989-2004 (", var.1[1], "%)", sep=""), cex=0.8)

z  <- eig.2.1
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim.2[1], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 1989-2004 (", var.1[2], "%)", sep=""), cex=0.8)

#############
# Era 2
z  <- eig.1.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim.1[2], lim.1[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF1 2005-2024 (", var.2[1], "%)", sep=""), cex=0.8)

z  <- eig.2.2
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim.2[1], lim.2[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext(paste("EOF2 2005-2024 (", var.2[2], "%)", sep=""), cex=0.8)


## so the loadings are very similar across the full time series and the two sub-eras

pc1_all <- as.matrix(SLP.NDJFM) %*% EOF.all$U[,1]
pc2_all <- as.matrix(SLP.NDJFM) %*% EOF.all$U[,2]


# average by winter year
m.NDJFM <- m[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec")]


win.yr <- if_else(m.NDJFM %in% c("Nov", "Dec"), as.numeric(as.character(yr.NDJFM))+1, as.numeric(as.character(yr.NDJFM))) 

pc_all_time_series <- data.frame(winter_year = 1948:2025,
                                 PC1 = tapply(pc1_all, win.yr, mean),
                                 PC2 = tapply(pc2_all, win.yr, mean))


plot_pc <- pc_all_time_series %>%
  pivot_longer(cols = -winter_year)

ggplot(plot_pc, aes(winter_year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(4,2)])

# out of curiosity - rolling 15 year pc1-pc2 correlations

pc_all_time_series$correlation <- pc_all_time_series$PC1_sd <- pc_all_time_series$PC2_sd <- 
  pc_all_time_series$PC1_AR1 <- pc_all_time_series$PC2_AR1 <-NA

for(i in 8:(nrow(pc_all_time_series)-7)){
  # i <- 8
  window <- (i-7):(i+7)
  
  pc_all_time_series$correlation[i] <- 
    cor(pc_all_time_series$PC1[window], pc_all_time_series$PC2[window])
  
  pc_all_time_series$PC1_sd[i] <- 
    sd(pc_all_time_series$PC1[window])
  
  pc_all_time_series$PC2_sd[i] <- 
    sd(pc_all_time_series$PC2[window])
  
  pc_all_time_series$PC1_AR1[i] <- 
    acf(pc_all_time_series$PC1[window], lag.max = 1, plot = FALSE)$acf[2]
  
  pc_all_time_series$PC2_AR1[i] <- 
    acf(pc_all_time_series$PC2[window], lag.max = 1, plot = FALSE)$acf[2]
   
}

ggplot(pc_all_time_series, aes(winter_year, correlation)) +
  geom_line() +
  geom_hline(yintercept = 0, lty = 2)
  
ggsave("./Figures/rolling_15_yr_correlation_N_Pac_slp_pc1-2.png", width = 6, height = 4, units = 'in')

# plot AR
plot_ar <- pc_all_time_series %>%
  select(winter_year, PC1_AR1, PC2_AR1) %>%
  pivot_longer(cols = -winter_year)

ggplot(plot_ar, aes(winter_year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(4,2)])

# plot SD
plot_sd <- pc_all_time_series %>%
  select(winter_year, PC1_sd, PC2_sd) %>%
  pivot_longer(cols = -winter_year)

ggplot(plot_sd, aes(winter_year, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(4,2)])

## now regressions

## plot regression of cellwise SST against mean SLP at center of EOF1 for SLP

# get SLP in EOF1 box
EOF1.x <- c(186.25, 208.75, 208.75, 186.25) 
EOF1.y <- c(56.25, 56.25, 46, 46)

xp <- cbind(EOF1.x, EOF1.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

SLP.box.NDJFM <- SLP.NDJFM
SLP.box.NDJFM[,!check] <- NA

SLP.box.monthly.mean <- rowMeans(SLP.box.NDJFM, na.rm = T)
SLP.box.winter.mean <- as.vector(scale(tapply(SLP.box.monthly.mean, win.yr, mean)))
names(SLP.box.winter.mean) <- 1948:2025
SLP.box.winter.mean <- SLP.box.winter.mean[names(SLP.box.winter.mean) %in% 1948:2024]

## load ERSST for the N. Pacific
# https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_31a3_72d5_401e.nc?sst[(1948-01-15):1:(2024-12-15)][(20):1:(68)][(130):1:(250)]

nc <- nc_open("./data/hawaii_soest_31a3_72d5_401e_d0ad_fd6c_67ce.nc")

# extract dates
ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

# extract study area
# 20-70 deg. N, 120-250 deg. E
x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

x; y


SST <- ncvar_get(nc, "sst", verbose = F)

# Change data from a 3-D array to a matrix of monthly data by grid point:
# First, reverse order of dimensions ("transpose" array)
SST <- aperm(SST, 3:1)  

# Change to matrix with column for each grid point, rows for monthly means
SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))  

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y))   
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


m <- months(d)  # Extracts months from the date vector
yr <- years(d)

# and plot 
SST.mean <- colMeans(SST)
z <- t(matrix(SST.mean,length(y)))  # Re-shape to a matrix with latitudes in columns, longitudes in rows
image(x,y,z, col=new.col)
contour(x, y, z, add=T, col="white")  
map('world2Hires',fill=F,add=T, lwd=2)

# identify columns in SST matrix corresponding to land
land <- is.na(colMeans(SST)) 

# For analysis, we only use the columns of the matrix with non-missing values:
X <- SST[,!land]

# remove seasonal means
f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(d)/12),]  # replicate means matrix for each year at each location

mu <- mu[rep(1:12, floor(length(d)/12)),] 


anom <- X - mu   # compute matrix of anomalies



# now detrend
anom.detr <- anom
for(i in 1:ncol(anom)) {
  # i <- 1
  xx = seq(1,nrow(anom))
  anom.detr[,i] = anom[,i] - predict(lm(anom[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}


# restrict to NDJFM
SST.NDJFM <- anom.detr[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec"),]

## get winter means for each cell -------------
SST.winter <- matrix(NA, ncol = ncol(SST.NDJFM), nrow = length(1948:2025))
colnames(SST.winter) <- colnames(SST.NDJFM)
rownames(SST.winter) <- 1948:2025

for(j in 1:ncol(SST.NDJFM)){
  # j <- 1
  
  SST.winter[,j] <- tapply(SST.NDJFM[,j], win.yr, mean) 
  
}

# remove 2025
SST.winter <- SST.winter[rownames(SST.winter) %in% 1948:2024,]

# now detrend
SST.winter.detr <- SST.winter

for(i in 1:ncol(SST.winter)) {
  # i <- 1
  xx = seq(1,nrow(SST.winter))
  SST.winter.detr[,i] = SST.winter[,i] - predict(lm(SST.winter[,i]~as.numeric(xx), na.action="na.exclude"), newdata=data.frame(xx=xx))
}


# regress 

regr.all <- regr.89.04 <- regr.05.24 <- NA

for(j in 1:ncol(SST.winter.detr)){
  # j <- 1
  
  mod <- lm(SST.winter.detr[,j] ~ SLP.box.winter.mean)
  regr.all[j] <- mod$coefficients[2]
  
  mod <- lm(SST.winter.detr[rownames(SST.winter.detr) %in% 1989:2004,j] ~ SLP.box.winter.mean[names(SLP.box.winter.mean) %in% 1989:2004])
  regr.89.04[j] <- mod$coefficients[2]  
  
  mod <- lm(SST.winter.detr[rownames(SST.winter.detr) %in% 2005:2024,j] ~ SLP.box.winter.mean[names(SLP.box.winter.mean) %in% 2005:2024])
  regr.05.24[j] <- mod$coefficients[2]   
  
}



# plot
lim <- range(regr.all, regr.89.04, regr.05.24)

png("./Figures/SST_fields_vs_SLP_box.png", width = 6, height = 12, units = 'in', res = 300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(3,1), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# regr all
z <- colMeans(SST)
z[!land]  <- regr.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
polygon(x = c(186.25, 208.75, 208.75, 186.25), y = c(56.25, 56.25, 46, 46), border = "white", lwd = 2)
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SST ~ SLP, 1948-2024", cex=0.8)

# regr 89-04
z <- colMeans(SST)
z[!land] <- regr.89.04
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SST ~ SLP, 1989-2004", cex=0.8)

# regr 05-24
z <- colMeans(SST)
z[!land] <- regr.05.24
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SST ~ SLP, 2005-2024", cex=0.8)

dev.off()

## only full time series
lim <- range(regr.all)
lim
png("./Figures/SST_fields_vs_SLP_box_full_time_series.png", width = 7, height = 6, units = 'in', res = 300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# regr all
z <- colMeans(SST)
z[!land]  <- regr.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")
polygon(x = c(186.25, 208.75, 208.75, 186.25), y = c(56.25, 56.25, 46, 46), border = "white", lwd = 2)
# mtext("Winter SST regressed on Aleutian Low", cex=0.8)

dev.off()

################################
# load goa sst
goa_sst <- read.csv(paste0(dir, "Data/goa.monthlySSTanomalies.csv"), row.names = 1) %>%
  filter(Month %in% c(1:3, 11, 12),
         Year %in% 1948:2024) %>%
  mutate(win.yr = if_else(Month %in% 11:12, Year+1, Year))

goa_win_sst <- data.frame(year = 1948:2024,
                          mean_anom = tapply(goa_sst$month.anom, goa_sst$win.yr, mean))
# detrend and scale

goa_win_sst$detr_anom <-
  goa_win_sst$mean_anom - predict(lm(goa_win_sst$mean_anom ~ goa_win_sst$year))



goa_win_sst$sc_detr_anom <- scale(as.vector(goa_win_sst$detr_anom))


# load ebs sst
ebs_sst <- read.csv("./data/ebs.monthlySSTanomalies.csv", row.names = 1) %>%
  filter(Month %in% c(1:3, 11, 12),
         Year %in% 1948:2024) %>%
  mutate(win.yr = if_else(Month %in% 11:12, Year+1, Year))

ebs_win_sst <- data.frame(year = 1948:2024,
                          mean_anom = tapply(ebs_sst$month.anom, ebs_sst$win.yr, mean))

# detrend and scale
ebs_win_sst$detr_anom <-
  ebs_win_sst$mean_anom - predict(lm(ebs_win_sst$mean_anom ~ ebs_win_sst$year))



ebs_win_sst$sc_detr_anom <- scale(as.vector(ebs_win_sst$detr_anom))

SLP.winter <- scale(SLP.winter)

# regress 

goa.all <- goa.89.04 <- goa.05.24 <- ebs.all <- ebs.89.04 <- ebs.05.24 <- NA

for(j in 1:ncol(SLP.winter)){
  # j <- 1
  
  # goa
  mod <- lm(goa_win_sst$sc_detr_anom ~ SLP.winter[rownames(SLP.winter) %in% 1948:2024,j])
  goa.all[j] <- mod$coefficients[2]
  
  mod <- lm(goa_win_sst$sc_detr_anom[goa_win_sst$year %in% 1989:2004]
            ~ SLP.winter[rownames(SLP.winter) %in% 1989:2004,j])
  goa.89.04[j] <- mod$coefficients[2]  

  mod <- lm(goa_win_sst$sc_detr_anom[goa_win_sst$year %in% 2005:2024]
            ~ SLP.winter[rownames(SLP.winter) %in% 2005:2024,j])
  goa.05.24[j] <- mod$coefficients[2]   
  
  # ebs
  mod <- lm(ebs_win_sst$sc_detr_anom ~ SLP.winter[rownames(SLP.winter) %in% 1948:2024,j])
  ebs.all[j] <- mod$coefficients[2]
  
  mod <- lm(ebs_win_sst$sc_detr_anom[ebs_win_sst$year %in% 1989:2004]
            ~ SLP.winter[rownames(SLP.winter) %in% 1989:2004,j])
  ebs.89.04[j] <- mod$coefficients[2]  
  
  mod <- lm(ebs_win_sst$sc_detr_anom[ebs_win_sst$year %in% 2005:2024]
            ~ SLP.winter[rownames(SLP.winter) %in% 2005:2024,j])
  ebs.05.24[j] <- mod$coefficients[2]  
  
}

# plot
lim <- range(goa.all, goa.89.04, goa.05.24, ebs.all, ebs.89.04, ebs.05.24)
lim <- -1
# goa all
z  <- goa.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 1948-2024", cex=0.8)

# goa 89-04
z  <- goa.89.04
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 1989-2004", cex=0.8)

# goa 05-24
z  <- goa.05.24
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 2005-2024", cex=0.8)

# ebs all
z  <- ebs.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1948-2024", cex=0.8)

# ebs 89-04
z  <- ebs.89.04
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1989-2004", cex=0.8)

# ebs 05-24
z  <- ebs.05.24
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 2005-2024", cex=0.8)

### 
# reverse regression

goa.all.r <- goa.89.04.r <- goa.05.24.r <- ebs.all.r <- ebs.89.04.r <- ebs.05.24.r <- NA

for(j in 1:ncol(SLP.winter)){
  # j <- 1
  
  # goa
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 1948:2024,j] ~ goa_win_sst$sc_detr_anom)
  goa.all.r[j] <- mod$coefficients[2]
  
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 1989:2004,j] ~
            goa_win_sst$sc_detr_anom[goa_win_sst$year %in% 1989:2004])
  goa.89.04.r[j] <- mod$coefficients[2]  
  
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 2005:2024,j] ~
              goa_win_sst$sc_detr_anom[goa_win_sst$year %in% 2005:2024])
  goa.05.24.r[j] <- mod$coefficients[2]   
  
  # ebs
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 1948:2024,j] ~ 
              ebs_win_sst$sc_detr_anom)
  ebs.all.r[j] <- mod$coefficients[2]
  
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 1989:2004,j] ~ 
              ebs_win_sst$sc_detr_anom[ebs_win_sst$year %in% 1989:2004])
  ebs.89.04.r[j] <- mod$coefficients[2]  
  
  mod <- lm(SLP.winter[rownames(SLP.winter) %in% 2005:2024,j] ~ 
              ebs_win_sst$sc_detr_anom[ebs_win_sst$year %in% 2005:2024])
  ebs.05.24.r[j] <- mod$coefficients[2]  
  
}

# plot
lim <- range(goa.all.r, goa.89.04.r, goa.05.24.r, ebs.all.r, ebs.89.04.r, ebs.05.24.r)

# goa all
z  <- goa.all.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 1948-2024", cex=0.8)

# goa 89-04
z  <- goa.89.04.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 1989-2004", cex=0.8)

# goa 05-24
z  <- goa.05.24.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 2005-2024", cex=0.8)

# ebs all
z  <- ebs.all.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1948-2024", cex=0.8)

# ebs 89-04
z  <- ebs.89.04.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1989-2004", cex=0.8)

# ebs 05-24
z  <- ebs.05.24.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 2005-2024", cex=0.8)

##
