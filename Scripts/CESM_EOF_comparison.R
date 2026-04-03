# PURPOSE: Compare ERA5 and CESM2 FCM/MDM EOF1 loadings for SST and SLP
# Author: Mike Litzow
#
# Produces four 4x3 panel plots (12 panels each):
#   Panel 1     = ERA5 observations
#   Panels 2-12 = 11 CESM2 ensemble members
#
# Plot 1: SST EOF1 — ERA5 obs vs FCM members
# Plot 2: SST EOF1 — ERA5 obs vs MDM members
# Plot 3: SLP EOF1 — ERA5 obs vs FCM members
# Plot 4: SLP EOF1 — ERA5 obs vs MDM members
#
# EOF fitting details:
#   - All months (not seasonal subset)
#   - Monthly anomalies from 1950-1979 climatology
#   - Linear detrending of each grid cell
#   - Latitude weighting: sqrt(cos(lat))
#   - Fitted via irlba (memory-efficient; avoids forming cells x cells matrix)
#   - Only leading EOF extracted
#   - Common period: 1950-2014

source("./Scripts/load.libs.functions.R")
library(irlba)

# ---- SETTINGS ----
yr.range   <- c(1950, 2014)
clim.range <- c(1950, 1979)
n.members  <- 11
lat.min    <- 20;  lat.max    <- 66
lon.min360 <- 110; lon.max360 <- 250

mapWorld <- map_data("world", wrap = c(20, 380))

# ---- HELPER FUNCTIONS ----

# Compute EOF1 from a raw monthly [lon x lat x time] array.
# Steps: trim to yr.range, anomaly from clim.range climatology,
#        detrend each cell, latitude-weight, irlba for leading EOF.
# Returns data frame: lon, lat, loading.
compute_eof1 <- function(arr, lons, lats, years, months) {

  # Trim to analysis period
  tidx   <- which(years >= yr.range[1] & years <= yr.range[2])
  arr    <- arr[, , tidx]
  yrs    <- years[tidx]
  mons   <- months[tidx]

  nlon <- length(lons); nlat <- length(lats); nt <- dim(arr)[3]

  # Monthly climatology from clim.range
  clim <- array(NA_real_, dim = c(nlon, nlat, 12))
  for (m in 1:12) {
    ci <- which(mons == m & yrs >= clim.range[1] & yrs <= clim.range[2])
    if (length(ci) > 0)
      clim[, , m] <- apply(arr[, , ci, drop = FALSE], c(1, 2),
                           mean, na.rm = TRUE)
  }

  # Anomalies
  anom <- arr
  for (t in seq_len(nt))
    anom[, , t] <- arr[, , t] - clim[, , mons[t]]

  # Reshape to [time x cells]
  mat  <- matrix(anom, nrow = nt)

  # Retain cells with data in > 50% of time steps
  good <- which(colSums(!is.na(mat)) > nt * 0.5)
  mat2 <- mat[, good]

  # Detrend each cell (linear trend on anomaly time series)
  mat2 <- apply(mat2, 2, function(x) {
    ok <- !is.na(x)
    if (sum(ok) > 5)
      x[ok] <- residuals(lm(x[ok] ~ seq_len(sum(ok))))
    x
  })

  # Latitude weights: sqrt(cos(lat))
  grid    <- expand.grid(lon = lons, lat = lats)
  weights <- sqrt(cos(grid$lat[good] * pi / 180))

  # Weight columns of data matrix — equivalent to weighted covariance SVD
  # but avoids forming the [cells x cells] covariance matrix
  mat.w <- sweep(mat2, 2, weights, "*")
  sv    <- irlba(mat.w, nv = 1, nu = 0)

  loading       <- rep(NA_real_, nlon * nlat)
  loading[good] <- sv$v[, 1]

  df          <- expand.grid(lon = lons, lat = lats)
  df$loading  <- loading
  df
}

# Sign convention: flip so dominant lobe is negative
orient_eof <- function(df) {
  if (mean(df$loading, na.rm = TRUE) > 0) df$loading <- -df$loading
  df
}

to360 <- function(lon) ifelse(lon < 0, lon + 360, lon)

# Single map panel
eof_panel <- function(df, title, col.lim) {
  df$lon360 <- to360(df$lon)
  ggplot() +
    geom_tile(data = filter(df, !is.na(loading)),
              aes(x = lon360, y = lat, fill = loading)) +
    geom_polygon(data = mapWorld,
                 aes(x = long, y = lat, group = group),
                 fill = "gray80", color = "gray50", linewidth = 0.15) +
    scale_fill_gradientn(
      colors   = rev(oceColorsPalette(64)),
      limits   = c(-col.lim, col.lim),
      na.value = "white",
      name     = "EOF1\nLoading"
    ) +
    coord_cartesian(xlim = c(110, 250), ylim = c(20, 66)) +
    labs(title = title) +
    theme_bw(base_size = 8) +
    theme(
      plot.title        = element_text(size = 7, hjust = 0.5),
      axis.title        = element_blank(),
      axis.text         = element_text(size = 5),
      panel.grid        = element_blank(),
      legend.position   = "right",
      legend.key.height = unit(0.8, "cm")
    )
}

# Build 4x3 panel plot: obs + 11 model members
make_12panel <- function(obs.df, member.dfs, var.label, model.label) {
  all.load <- c(obs.df$loading, unlist(lapply(member.dfs, `[[`, "loading")))
  col.lim  <- max(abs(all.load), na.rm = TRUE)

  panels <- c(
    list(eof_panel(obs.df, paste0("ERA5 Obs (", var.label, ")"), col.lim)),
    lapply(seq_along(member.dfs), function(i)
      eof_panel(member.dfs[[i]], paste0(model.label, " member ", i), col.lim))
  )

  wrap_plots(panels, ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0(
        "EOF1 Loadings — ERA5 vs CESM2 ", model.label, " (", var.label, ")\n",
        "All months, 1950-2014, anomaly from 1950-1979 climatology"
      ),
      theme = theme(plot.title = element_text(size = 10, hjust = 0.5))
    )
}

# Load one CESM2 member; returns list with arr, lons, lats, years, months
load_cesm <- function(file, varname) {
  nc  <- nc_open(file)
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  dat <- ncvar_get(nc, varname)
  nc_close(nc)
  if (length(dim(dat)) == 4) dat <- dat[, , 1, ]   # drop z_t if present
  dat[dat > 1e20] <- NA

  lat.idx <- which(lat >= lat.min & lat <= lat.max)
  lon.idx <- which(lon >= lon.min360 & lon <= lon.max360)
  dat     <- dat[lon.idx, lat.idx, ]

  nt    <- dim(dat)[3]
  dates <- seq(as.Date("1850-01-15"), by = "month", length.out = nt)
  list(arr    = dat,
       lons   = lon[lon.idx],
       lats   = lat[lat.idx],
       years  = as.integer(format(dates, "%Y")),
       months = as.integer(format(dates, "%m")))
}

# ---- ERA5 SST ----
message("Loading ERA5 SST...")
nc       <- nc_open("./Data/era5_sst_NP_1950_2025.nc")
sst.raw  <- ncvar_get(nc, "sst") - 273.15
lons.sst <- ncvar_get(nc, "longitude")
lats.sst <- ncvar_get(nc, "latitude")
time     <- ncvar_get(nc, "valid_time")
nc_close(nc)

dates.e <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
yrs.e   <- as.integer(format(dates.e, "%Y"))
mons.e  <- as.integer(format(dates.e, "%m"))

message("Computing ERA5 SST EOF1...")
eof.sst.obs <- orient_eof(
  compute_eof1(sst.raw, lons.sst, lats.sst, yrs.e, mons.e)
)
rm(sst.raw); gc()

# ---- ERA5 SLP ----
message("Loading ERA5 SLP...")
nc       <- nc_open("./Data/era5_slp_NP_1950_2025.nc")
slp.raw  <- ncvar_get(nc, "msl") / 100     # Pa -> hPa
lons.slp <- ncvar_get(nc, "longitude")
lats.slp <- ncvar_get(nc, "latitude")
time     <- ncvar_get(nc, "valid_time")
nc_close(nc)

dates.s <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
yrs.s   <- as.integer(format(dates.s, "%Y"))
mons.s  <- as.integer(format(dates.s, "%m"))

# Restrict to North Pacific domain
lat.idx.slp <- which(lats.slp >= lat.min & lats.slp <= lat.max)
lon.idx.slp <- which(to360(lons.slp) >= lon.min360 &
                       to360(lons.slp) <= lon.max360)
slp.crop  <- slp.raw[lon.idx.slp, lat.idx.slp, ]
lons.slp2 <- lons.slp[lon.idx.slp]
lats.slp2 <- lats.slp[lat.idx.slp]
rm(slp.raw); gc()

message("Computing ERA5 SLP EOF1...")
eof.slp.obs <- orient_eof(
  compute_eof1(slp.crop, lons.slp2, lats.slp2, yrs.s, mons.s)
)
rm(slp.crop); gc()

# ---- CESM2 FILE LISTS ----
fcm.sst.files <- head(sort(list.files("./Data/CESM2 ensemble/SST/FCM",
                                      pattern = "\\.nc$", full.names = TRUE)),
                      n.members)
mdm.sst.files <- head(sort(list.files("./Data/CESM2 ensemble/SST/MDM",
                                      pattern = "\\.nc$", full.names = TRUE)),
                      n.members)
fcm.slp.files <- head(sort(list.files("./Data/CESM2 ensemble/SLP/FCM",
                                      pattern = "\\.nc$", full.names = TRUE)),
                      n.members)
mdm.slp.files <- head(sort(list.files("./Data/CESM2 ensemble/SLP/MDM",
                                      pattern = "\\.nc$", full.names = TRUE)),
                      n.members)

# Process one CESM2 member: load, compute EOF1, orient
eof_cesm <- function(file, varname, scale.pa = FALSE) {
  d <- load_cesm(file, varname)
  if (scale.pa) d$arr <- d$arr / 100
  orient_eof(compute_eof1(d$arr, d$lons, d$lats, d$years, d$months))
}

# ---- PROCESS CESM2 MEMBERS ----
message("Processing FCM SST members (", length(fcm.sst.files), ")...")
fcm.sst.eofs <- lapply(seq_along(fcm.sst.files), function(i) {
  message("  member ", i); eof_cesm(fcm.sst.files[i], "SST")
})

message("Processing MDM SST members (", length(mdm.sst.files), ")...")
mdm.sst.eofs <- lapply(seq_along(mdm.sst.files), function(i) {
  message("  member ", i); eof_cesm(mdm.sst.files[i], "SST")
})

message("Processing FCM SLP members (", length(fcm.slp.files), ")...")
fcm.slp.eofs <- lapply(seq_along(fcm.slp.files), function(i) {
  message("  member ", i); eof_cesm(fcm.slp.files[i], "PSL", scale.pa = TRUE)
})

message("Processing MDM SLP members (", length(mdm.slp.files), ")...")
mdm.slp.eofs <- lapply(seq_along(mdm.slp.files), function(i) {
  message("  member ", i); eof_cesm(mdm.slp.files[i], "PSL", scale.pa = TRUE)
})

# ---- PLOTS ----
message("Plotting...")

p1 <- make_12panel(eof.sst.obs, fcm.sst.eofs, "SST", "FCM")
p2 <- make_12panel(eof.sst.obs, mdm.sst.eofs, "SST", "MDM")
p3 <- make_12panel(eof.slp.obs, fcm.slp.eofs, "SLP", "FCM")
p4 <- make_12panel(eof.slp.obs, mdm.slp.eofs, "SLP", "MDM")

print(p1)
print(p2)
print(p3)
print(p4)
