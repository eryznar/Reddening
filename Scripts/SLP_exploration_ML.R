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

SLP.raw <- SLP[m %in% c("Jan", "Feb", "Mar", "Nov", "Dec"),]

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

## get winter means for each cell -------------
SLP.winter <- matrix(NA, ncol = ncol(SLP.NDJFM), nrow = length(1948:2025))
colnames(SLP.winter) <- colnames(SLP.NDJFM)
rownames(SLP.winter) <- 1948:2025

for(j in 1:ncol(SLP.NDJFM)){
  # j <- 1
  
  SLP.winter[,j] <- tapply(SLP.NDJFM[,j], win.yr, mean) 
  
}

View(SLP.winter)

## winter raw SLP

# SLP.winter.raw <- matrix(NA, ncol = ncol(SLP.raw), nrow = length(1948:2025))
# colnames(SLP.winter.raw) <- colnames(SLP.raw)
# rownames(SLP.winter.raw) <- 1948:2025
# 
# for(j in 1:ncol(SLP.raw)){
#   # j <- 1
#   
#   SLP.winter.raw[,j] <- tapply(SLP.raw[,j], win.yr, mean) 
#   
# }
# 
# View(SLP.winter.raw)
# 
# SLP.winter.raw <- scale(SLP.winter.raw)

# load goa sst
goa_sst <- read.csv("./data/goa.monthlySSTanomalies.csv", row.names = 1) %>%
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

# goa all
z  <- goa.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 1948-2024", cex=0.8)

# goa 89-04
z  <- goa.89.04
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 1989-2004", cex=0.8)

# goa 05-24
z  <- goa.05.24
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("GOA SST ~ SLP, 2005-2024", cex=0.8)

# ebs all
z  <- ebs.all
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1948-2024", cex=0.8)

# ebs 89-04
z  <- ebs.89.04
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1989-2004", cex=0.8)

# ebs 05-24
z  <- ebs.05.24
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(lim[1], -lim[1]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
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
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 1948-2024", cex=0.8)

# goa 89-04
z  <- goa.89.04.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 1989-2004", cex=0.8)

# goa 05-24
z  <- goa.05.24.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("SLP ~ GOA SST, 2005-2024", cex=0.8)

# ebs all
z  <- ebs.all.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1948-2024", cex=0.8)

# ebs 89-04
z  <- ebs.89.04.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 1989-2004", cex=0.8)

# ebs 05-24
z  <- ebs.05.24.r
z <- t(matrix(z, length(y))) 
image(x,y,z, col=new.col, xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=c(-lim[2], lim[2]))#, legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))
contour(x, y, z, add=T, drawlabels = F, lwd=0.7, col="grey") 
map('world2Hires', fill=F, add=T, lwd=1)
mtext("ebs SST ~ SLP, 2005-2024", cex=0.8)

##
