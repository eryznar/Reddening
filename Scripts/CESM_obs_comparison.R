# PURPOSE: Compare observed and modeled AR1 values for SST in the GOA and EBS
# Author: Mike Litzow

# LOAD LIBS/FUNCTIONS ----------------------------------
source("./Scripts/load.libs.functions.R")

# MAP OF GOA AND EBS POLYGONS ----------------------------------

# Pacific-centered world map (longitudes wrapped to 20-380 range)
mapWorld <- map_data("world", wrap = c(20, 380))

# Define GOA polygon in 0-360 space (for masking, convert to -180/180 below)
goa.x <- c(198, 198, 203, 205, 208, 225, 231, 201)
goa.y <- c(54, 55.5, 57.5, 59, 61, 61, 54, 54)

goa.poly <- data.frame(lon = goa.x, lat = goa.y, region = "GOA")

# Define EBS polygon in 0-360 space
ebs.x <- c(183, 183, 203, 203, 191)
ebs.y <- c(53, 65, 65, 57.5, 53)

ebs.poly <- data.frame(lon = ebs.x, lat = ebs.y, region = "EBS")

# Convert to -180/180 for sf masking
goa.x <- ifelse(goa.x > 180, goa.x - 360, goa.x)
ebs.x <- ifelse(ebs.x > 180, ebs.x - 360, ebs.x)

# Combine GOA and EBS polygons (0-360 lons for Pacific-centered plot)
polys <- bind_rows(goa.poly, ebs.poly)

# North Pacific ERA5 domain (20-66N, 110-250E) — single rectangle in 0-360 space
np.rect <- data.frame(xmin = 110, xmax = 250, ymin = 20, ymax = 66)

# Plot
map.plot <- ggplot() +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_rect(data = np.rect,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "forestgreen", color = "forestgreen", alpha = 0.1, linewidth = 0.7) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region, fill = region, color = region),
               alpha = 0.4, linewidth = 1) +
  scale_fill_manual(values = c("GOA" = "steelblue", "EBS" = "darkorange")) +
  scale_color_manual(values = c("GOA" = "steelblue4", "EBS" = "darkorange4")) +
  coord_cartesian(xlim = c(100, 260), ylim = c(15, 70)) +
  labs(title = "North Pacific ERA5 Domain with GOA and EBS Study Regions",
       x = "Longitude", y = "Latitude",
       fill = "Region", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

print(map.plot)

# DOWNLOAD ERA5 SST ----------------------------------

library(ecmwfr)

# Set CDS API key (run once, then comment out — stores credentials in macOS Keychain)
# Find your UID and API key at: https://cds.climate.copernicus.eu -> your profile
# wf_set_key(user = "YOUR_CDS_UID", key = "YOUR_CDS_API_KEY", service = "cds")

# Full North Pacific domain: 20-66N, 110-250E (area: N, W, S, E in -180/180)
request <- list(
  dataset_short_name = "reanalysis-era5-single-levels-monthly-means",
  product_type    = "monthly_averaged_reanalysis",
  variable        = "sea_surface_temperature",
  year            = as.character(1950:2024),
  month           = sprintf("%02d", 1:12),
  time            = "00:00",
  area            = c(66, 110, 20, -110),    # N, W, S, E (250E = -110W)
  data_format     = "netcdf",
  download_format = "unarchived",
  target          = "era5_sst_NP_1950_2024.nc"
)

# Slow down polling to avoid CDS rate limit
options(ecmwfr.sleep = 120)   # wait 120s between status checks

wf_request(
  request  = request,
  transfer = TRUE,
  path     = "./Data"
)

# CALCULATE OBSERVED AR1 AND SD FROM ERA5 ----------------------------------

# Load ERA5 SST netCDF (North Pacific domain)
nc <- nc_open("./Data/era5_sst_NP_1950_2024.nc")

sst  <- ncvar_get(nc, "sst")        # [longitude x latitude x valid_time], units: K
lons <- ncvar_get(nc, "longitude")
lats <- ncvar_get(nc, "latitude")
time <- ncvar_get(nc, "valid_time") # seconds since 1970-01-01
nc_close(nc)

# Convert time to dates
dates  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
months <- as.integer(format(dates, "%m"))
years  <- as.integer(format(dates, "%Y"))

# Convert SST from Kelvin to Celsius
sst <- sst - 273.15

# Build polygon sf objects for masking (close rings by repeating first coordinate)
goa.coords <- rbind(cbind(goa.x, goa.y), c(goa.x[1], goa.y[1]))
ebs.coords <- rbind(cbind(ebs.x, ebs.y), c(ebs.x[1], ebs.y[1]))

goa.sf <- st_sf(geometry = st_sfc(st_polygon(list(goa.coords)), crs = 4326))
ebs.sf <- st_sf(geometry = st_sfc(st_polygon(list(ebs.coords)), crs = 4326))

# Build grid of all lon/lat combinations matching array order
grid    <- expand.grid(lon = lons, lat = lats)
grid.sf <- st_as_sf(grid, coords = c("lon", "lat"), crs = 4326)

# Identify which cells fall inside each polygon / domain
in.goa <- lengths(st_within(grid.sf, goa.sf)) > 0
in.ebs <- lengths(st_within(grid.sf, ebs.sf)) > 0
in.np  <- !is.na(as.vector(sst[,,1]))   # all ocean cells

# Latitude-based weights (cos(lat) proportional to cell area)
weights <- cos(grid$lat * pi / 180)

# --- Step 1: area-wide weighted mean SST for each domain ---
weighted_region_mean <- function(mask) {
  w <- weights[mask]
  apply(sst, 3, function(slice) {
    vals <- as.vector(slice)[mask]
    if (all(is.na(vals))) return(NA_real_)
    weighted.mean(vals, w, na.rm = TRUE)
  })
}

goa.sst <- weighted_region_mean(in.goa)
ebs.sst <- weighted_region_mean(in.ebs)
np.sst  <- weighted_region_mean(in.np)

# Assemble monthly data frame
monthly <- data.frame(
  year  = years,
  month = months,
  GOA   = goa.sst,
  EBS   = ebs.sst,
  NP    = np.sst
)

# --- Step 2: monthly climatology for each domain (1950-1979) ---
clim <- monthly %>%
  filter(year >= 1950, year <= 1979) %>%
  group_by(month) %>%
  summarise(GOA.clim = mean(GOA, na.rm = TRUE),
            EBS.clim = mean(EBS, na.rm = TRUE),
            NP.clim  = mean(NP,  na.rm = TRUE),
            .groups = "drop")

# --- Step 3: anomaly time series (°C) for full period ---
monthly <- monthly %>%
  left_join(clim, by = "month") %>%
  mutate(GOA.anom = GOA - GOA.clim,
         EBS.anom = EBS - EBS.clim,
         NP.anom  = NP  - NP.clim)

# --- Step 4: detrend each domain's anomaly time series ---
detrend <- function(x) residuals(lm(x ~ seq_along(x)))

monthly <- monthly %>%
  mutate(GOA.anom = detrend(GOA.anom),
         EBS.anom = detrend(EBS.anom),
         NP.anom  = detrend(NP.anom))

# --- Step 5: restrict to winter (Nov-Mar), year = January year ---
monthly.win <- monthly %>%
  filter(month %in% c(11, 12, 1, 2, 3)) %>%
  mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year))

# Winter mean anomaly per domain per year
winter <- monthly.win %>%
  group_by(win.year) %>%
  summarise(GOA = mean(GOA.anom, na.rm = TRUE),
            EBS = mean(EBS.anom, na.rm = TRUE),
            NP  = mean(NP.anom,  na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = win.year) %>%
  arrange(year)

# --- Step 6: 15-year rolling AR(1) and SD (right-aligned) ---
roll_ar1 <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2])
}
roll_sd <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = sd)
}

winter <- winter %>%
  mutate(GOA.ar1 = roll_ar1(GOA),
         EBS.ar1 = roll_ar1(EBS),
         NP.ar1  = roll_ar1(NP),
         GOA.sd  = roll_sd(GOA),
         EBS.sd  = roll_sd(EBS),
         NP.sd   = roll_sd(NP))

# --- Step 7: plot ---
region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4", "NP" = "forestgreen")

ar1.plot <- winter %>%
  select(year, GOA.ar1, EBS.ar1, NP.ar1) %>%
  pivot_longer(-year, names_to = "region", values_to = "ar1") %>%
  mutate(region = sub(".ar1", "", region, fixed = TRUE)) %>%
  ggplot(aes(x = year, y = ar1, color = region)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = region.colors) +
  labs(title = "15-year Rolling AR(1) — Winter SST Anomaly",
       x = "Year", y = "AR(1)", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

sd.plot <- winter %>%
  select(year, GOA.sd, EBS.sd, NP.sd) %>%
  pivot_longer(-year, names_to = "region", values_to = "sd") %>%
  mutate(region = sub(".sd", "", region, fixed = TRUE)) %>%
  ggplot(aes(x = year, y = sd, color = region)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = region.colors) +
  labs(title = "15-year Rolling SD — Winter SST Anomaly",
       x = "Year", y = "SD (°C)", color = "Region") +
  theme_bw() +
  theme(legend.position = "bottom")

print(ar1.plot / sd.plot)

# AUTOCORRELATION BY ERA ----------------------------------

# Use full monthly detrended anomaly time series (all months)
era.breaks <- list(
  "1950-1988" = c(1950, 1988),
  "1989-2000" = c(1989, 2000),
  "2001-2025" = c(2001, 2025)
)
era.colors <- c("1950-1988" = "steelblue4", "1989-2000" = "darkorange3", "2001-2025" = "firebrick3")

max.lag <- 60

# Compute ACF for one time series and return as data frame
acf_df <- function(x, max.lag) {
  a <- acf(x, lag.max = max.lag, plot = FALSE, na.action = na.pass)
  data.frame(lag = as.integer(a$lag), acf = as.numeric(a$acf))
}

# Build ACF data frame for all domains and eras
acf.out <- data.frame()

for (era in names(era.breaks)) {
  yr.range <- era.breaks[[era]]
  dat <- monthly %>%
    filter(year >= yr.range[1], year <= yr.range[2]) %>%
    arrange(year, month)

  for (dom in c("GOA", "EBS", "NP")) {
    col <- paste0(dom, ".anom")
    a   <- acf_df(dat[[col]], max.lag)
    acf.out <- rbind(acf.out,
                     data.frame(era = era, domain = dom, lag = a$lag, acf = a$acf))
  }
}

# Plot — one panel per domain
acf.plot <- ggplot(acf.out, aes(x = lag, y = acf, color = era)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray50") +
  geom_line() +
  geom_point(size = 1.2) +
  facet_wrap(~ domain, ncol = 1) +
  scale_color_manual(values = era.colors) +
  scale_x_continuous(breaks = seq(0, max.lag, by = 12)) +
  labs(title = "Autocorrelation of Monthly SST Anomaly by Era",
       x = "Lag (months)", y = "Autocorrelation", color = "Era") +
  theme_bw() +
  theme(legend.position = "bottom")

print(acf.plot)

# WAVELET ANALYSIS ----------------------------------

library(WaveletComp)

# Full monthly detrended anomaly time series, sorted chronologically
monthly.sorted <- monthly %>% arrange(year, month)

# Decimal year for x-axis
monthly.sorted <- monthly.sorted %>%
  mutate(dec.year = year + (month - 0.5) / 12)

domains <- c("GOA", "EBS", "NP")

# Run wavelet analysis and extract power/significance into a data frame
wavelet_df <- function(x, dec.year, dom) {
  wt <- analyze.wavelet(
    my.data      = data.frame(x = x),
    my.series    = "x",
    loess.span   = 0,
    dt           = 1/12,
    dj           = 1/20,
    lowerPeriod  = 0.5,
    upperPeriod  = 32,
    make.pval    = TRUE,
    n.sim        = 100,
    verbose      = FALSE
  )

  # wt$axis.1 indexes the non-padded time steps within the padded matrix
  power.df <- expand.grid(
    period = wt$Period,
    time   = dec.year
  )
  power.df$power  <- as.vector(wt$Power[, wt$axis.1])
  power.df$sig    <- as.vector(wt$Power.pval[, wt$axis.1]) < 0.05
  power.df$domain <- dom

  # Subset COI to non-padded time steps
  coi.df <- data.frame(time = dec.year, coi = wt$coi.1[wt$axis.1], domain = dom)

  list(power = power.df, coi = coi.df)
}

wt.list <- lapply(domains, function(dom) {
  wavelet_df(monthly.sorted[[paste0(dom, ".anom")]], monthly.sorted$dec.year, dom)
})

power.all <- do.call(rbind, lapply(wt.list, `[[`, "power"))
coi.all   <- do.call(rbind, lapply(wt.list, `[[`, "coi"))

# Mark cells outside the COI for shading
power.all <- power.all %>%
  left_join(coi.all %>% select(time, domain, coi), by = c("time", "domain")) %>%
  mutate(outside.coi = period > coi)

# Plot
wavelet.plot <- ggplot(power.all, aes(x = time, y = period, fill = log2(power))) +
  geom_raster(interpolate = TRUE) +
  # shade outside-COI region with semi-transparent white overlay
  geom_raster(data = filter(power.all, outside.coi),
              aes(x = time, y = period), fill = "white", alpha = 0.6,
              inherit.aes = FALSE) +
  # significance contour (only within COI)
  geom_contour(data = filter(power.all, !outside.coi),
               aes(x = time, y = period, z = as.integer(sig)),
               inherit.aes = FALSE,
               breaks = 0.5, color = "black", linewidth = 0.4) +
  scale_fill_distiller(palette = "Spectral", direction = -1,
                       name = expression(log[2](power))) +
  scale_y_log10(breaks = c(0.5, 1, 2, 4, 8, 16, 32),
                labels = c(0.5, 1, 2, 4, 8, 16, 32)) +
  facet_wrap(~ domain, ncol = 1) +
  labs(title = "Wavelet Power Spectrum — Monthly SST Anomaly (1950–2025)",
       x = "Year", y = "Period (years)") +
  theme_bw() +
  theme(legend.position = "right")

print(wavelet.plot)
