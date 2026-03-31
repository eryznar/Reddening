# PURPOSE: Compare observed and modeled AR1 values for SST in the GOA and EBS
# Author: Mike Litzow

# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# MAP OF GOA AND EBS POLYGONS ----------------------------------

mapWorld <- map_data("world")

# Define GOA polygon (convert 0-360 longitudes to -180 to 180)
goa.x <- c(198, 198, 203, 205, 208, 225, 231, 201)
goa.x <- ifelse(goa.x > 180, goa.x - 360, goa.x)
goa.y <- c(54, 55.5, 57.5, 59, 61, 61, 54, 54)

goa.poly <- data.frame(lon = goa.x, lat = goa.y, region = "GOA")

# Define EBS polygon (convert 0-360 longitudes to -180 to 180)
ebs.x <- c(183, 183, 203, 203, 191)
ebs.x <- ifelse(ebs.x > 180, ebs.x - 360, ebs.x)
ebs.y <- c(53, 65, 65, 57.5, 53)

ebs.poly <- data.frame(lon = ebs.x, lat = ebs.y, region = "EBS")

# Combine polygons
polys <- bind_rows(goa.poly, ebs.poly)

# Plot
map.plot <- ggplot() +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region, fill = region, color = region),
               alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("GOA" = "steelblue", "EBS" = "darkorange")) +
  scale_color_manual(values = c("GOA" = "steelblue4", "EBS" = "darkorange4")) +
  coord_cartesian(xlim = c(-180, -125), ylim = c(45, 70)) +
  labs(title = "GOA and EBS Study Regions",
       x = "Longitude", y = "Latitude",
       fill = "Region", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

print(map.plot)

# DOWNLOAD ERA5 SST ----------------------------------

library(ecmwfr)

# Set CDS API key (only needs to be run once per machine)
# wf_set_key(user = "your_uid", key = "your_api_key", service = "cds")

# Spatial domain encompassing GOA and EBS regions (N, W, S, E)
request <- list(
  product_type    = "monthly_averaged_reanalysis",
  variable        = "sea_surface_temperature",
  year            = as.character(1950:2024),
  month           = sprintf("%02d", 1:12),
  time            = "00:00",
  area            = c(70, -180, 45, -125),   # N, W, S, E
  data_format     = "netcdf",
  download_format = "unarchived"
)

wf_request(
  request  = request,
  transfer = TRUE,
  path     = "./Data",
  filename = "era5_sst_GOA_EBS_1950_2024.nc",
  dataset  = "reanalysis-era5-single-levels-monthly-means"
)
