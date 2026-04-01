# PURPOSE: Download ERA5 monthly mean SLP for the North Pacific domain
# Author: Mike Litzow

library(ecmwfr)

# Set CDS API key (run once, then comment out — stores credentials in macOS Keychain)
# Find your personal access token at: https://cds.climate.copernicus.eu -> your profile
# wf_set_key(key = "YOUR_PERSONAL_ACCESS_TOKEN", service = "cds")

# North Pacific domain: 20-70N, 120-250E (area: N, W, S, E in -180/180)
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  product_type       = "monthly_averaged_reanalysis",
  variable           = "mean_sea_level_pressure",
  year               = as.character(1950:2025),
  month              = sprintf("%02d", 1:12),
  time               = "00:00",
  area               = c(70, 120, 20, -110),   # N, W, S, E (250E = -110W)
  data_format        = "netcdf",
  download_format    = "unarchived",
  target             = "era5_slp_NP_1950_2025.nc"
)

# Slow down polling to avoid CDS rate limit
options(ecmwfr.sleep = 120)

wf_request(
  request  = request,
  transfer = TRUE,
  path     = "./Data"
)

# MAP OF WINTER MEAN SLP ----------------------------------

library(tidyverse)
library(ncdf4)

# Load SLP netCDF
nc  <- nc_open("./Data/era5_slp_NP_1950_2025.nc")
slp  <- ncvar_get(nc, "msl")         # [longitude x latitude x valid_time], units: Pa
lons <- ncvar_get(nc, "longitude")
lats <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "valid_time")  # seconds since 1970-01-01
nc_close(nc)

# Convert to hPa and extract dates
slp    <- slp / 100
dates  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
months <- as.integer(format(dates, "%m"))

# Winter (Nov-Mar) mean at each grid cell
win.idx <- which(months %in% c(11, 12, 1, 2, 3))
slp.win <- apply(slp[, , win.idx], c(1, 2), mean, na.rm = TRUE)

# Build data frame — convert lons > 180 to 0-360 for Pacific-centered plot
slp.df <- expand.grid(lon = lons, lat = lats)
slp.df$slp <- as.vector(slp.win)
slp.df <- slp.df %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon)) %>%
  filter(!is.na(slp))

# Pacific-centered basemap
mapWorld <- map_data("world", wrap = c(20, 380))

slp.map <- ggplot() +
  geom_raster(data = slp.df, aes(x = lon, y = lat, fill = slp)) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdBu", direction = -1,
                       name = "SLP (hPa)") +
  coord_cartesian(xlim = c(120, 250), ylim = c(20, 70)) +
  labs(title = "Winter Mean SLP (Nov–Mar, 1950–2025)",
       x = "Longitude", y = "Latitude") +
  theme_bw()

print(slp.map)

# EOF ANALYSIS OF SLP ----------------------------------

library(irlba)  # for fast truncated SVD (only compute 2 leading modes)

years  <- as.integer(format(dates, "%Y"))

# --- Step 1: monthly anomalies (subtract monthly climatology, 1950-1979) ---
n.lon  <- length(lons)
n.lat  <- length(lats)
n.time <- dim(slp)[3]

slp.anom <- array(NA_real_, dim = dim(slp))
base     <- years >= 1950 & years <= 1979

for (m in 1:12) {
  t.base <- which(base & months == m)
  t.all  <- which(months == m)
  clim   <- apply(slp[, , t.base, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  for (t in t.all) slp.anom[, , t] <- slp[, , t] - clim
}

# --- Step 2: detrend each cell's anomaly time series ---
X        <- cbind(1, seq_len(n.time))
H        <- X %*% solve(t(X) %*% X) %*% t(X)
IH       <- diag(n.time) - H
anom.mat <- matrix(slp.anom, nrow = n.lon * n.lat, ncol = n.time)  # [cells x time]

# All cells included (land + ocean); SLP is defined everywhere
detrend.mat <- anom.mat %*% t(IH)

# --- Step 3: area weighting by sqrt(cos(lat)) before covariance EOF ---
# Weight applied to each cell; sqrt because SVD of weighted matrix gives
# covariance-matrix eigenvectors when weights = sqrt(cell area)
grid <- expand.grid(lon = lons, lat = lats)
w    <- sqrt(cos(grid$lat * pi / 180))

weighted.mat <- sweep(detrend.mat, 1, w, "*")  # weight each cell's time series

# --- Step 4: truncated SVD — only compute 2 leading modes ---
# SVD of [cells x time]: U = spatial patterns, V = PC time series
sv <- irlba(weighted.mat, nv = 2)

# Recover unweighted EOF spatial patterns (divide back out area weight)
eof1 <- sv$u[, 1] / w
eof2 <- sv$u[, 2] / w

# PC time series (standardized)
pc1 <- sv$v[, 1]
pc2 <- sv$v[, 2]

# --- Step 5: spatial EOF plots ---
eof.df <- expand.grid(lon = lons, lat = lats) %>%
  mutate(lon  = ifelse(lon < 0, lon + 360, lon),
         EOF1 = eof1,
         EOF2 = eof2)

box1.df <- data.frame(x = c(191, 191, 208, 208, 191),
                      y = c(44,  55,  55,  44,  44))

eof1.map <- ggplot() +
  geom_raster(data = eof.df, aes(x = lon, y = lat, fill = EOF1)) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.2) +
  geom_path(data = box1.df, aes(x = x, y = y),
            color = "black", linewidth = 0.8) +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "Loading") +
  coord_cartesian(xlim = c(120, 250), ylim = c(20, 70)) +
  labs(title = "EOF1",
       x = "Longitude", y = "Latitude") +
  theme_bw()

eof2.map <- ggplot() +
  geom_raster(data = eof.df, aes(x = lon, y = lat, fill = EOF2)) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = NA, color = "gray30", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "Loading") +
  coord_cartesian(xlim = c(120, 250), ylim = c(20, 70)) +
  labs(title = "EOF2",
       x = "Longitude", y = "Latitude") +
  theme_bw()

print(eof1.map / eof2.map)

# --- Step 6: PC time series plots ---
pc.df <- data.frame(
  date = dates,
  year = years,
  month = months,
  PC1  = pc1,
  PC2  = pc2
)

pc.plot <- pc.df %>%
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(x = date, y = score)) +
  geom_line(linewidth = 0.4, color = "steelblue4") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  facet_wrap(~ PC, ncol = 1) +
  labs(title = "SLP PC Time Series (monthly, detrended anomalies)",
       x = "Year", y = "PC Score") +
  theme_bw()

print(pc.plot)

# ALEUTIAN LOW SLP TIME SERIES ----------------------------------

# Box coordinates (0-360 longitude space): 191-208E, 44-55N
# Identify grid cells inside the box
al.mask <- grid$lon >= 191 & grid$lon <= 208 &
           grid$lat >=  44 & grid$lat <=  55

# Mean detrended SLP anomaly inside box for each time step
al.slp <- apply(detrend.mat[al.mask, ], 2, mean, na.rm = TRUE)

# Assemble monthly data frame
al.monthly <- data.frame(
  year  = years,
  month = months,
  SLP   = al.slp
)

# Winter (Nov-Mar) means, year = January year
al.winter <- al.monthly %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
  group_by(win.year) %>%
  summarise(SLP = mean(SLP, na.rm = TRUE), .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year)

# Save to Output
write.csv(al.winter, "./Output/AL_winter_SLP_anomaly.csv", row.names = FALSE)
