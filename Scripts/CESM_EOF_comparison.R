# PURPOSE: Compare ERA5 and CESM2 FCM/MDM EOF1 loadings for SST and SLP
# Author: Mike Litzow
#
# Produces four 4x3 panel plots (12 panels each):
#   Panel 1     = ERA5 observations
#   Panels 2-12 = 11 CESM2 ensemble members
#
# EOF fitting details:
#   - All months (not seasonal subset)
#   - Monthly anomalies from 1950-1979 climatology
#   - Linear detrending of each grid cell
#   - Latitude weighting: sqrt(cos(lat))
#   - Temporal covariance approach: eigen([time x time]) — avoids memory issues
#     and is more numerically stable than SVD on large spatial matrices
#   - Only leading EOF extracted
#   - Common period: 1950-2014
#   - Same spatial domain for SST and SLP: 20-66N, 110-250E

source("./Scripts/load.libs.functions.R")

# ---- SETTINGS ----
yr.range   <- c(1950, 2014)
clim.range <- c(1950, 1979)
n.members  <- 11
lat.min    <- 20;  lat.max    <- 66
lon.min360 <- 110; lon.max360 <- 250

mapWorld <- map_data("world", wrap = c(20, 380))
to360    <- function(lon) ifelse(lon < 0, lon + 360, lon)

# ---- EOF COMPUTATION ----
# Input:  arr [lon x lat x time], coordinate vectors, year/month index vectors
# Method: temporal covariance eigen-decomposition
#   1. Trim to yr.range
#   2. Monthly anomaly relative to clim.range climatology
#   3. Linear detrend each cell
#   4. Latitude-weight columns: w = sqrt(cos(lat))
#   5. Form [time x time] covariance C = X_w X_w^T / (n-1)
#   6. Leading eigenvector of C = PC1 time series
#   7. Project X_w^T onto PC1 to recover EOF1 spatial loadings
# Returns data frame: lon, lat, loading

compute_eof1 <- function(arr, lons, lats, years, months) {

  tidx <- which(years >= yr.range[1] & years <= yr.range[2])
  arr  <- arr[, , tidx]
  yrs  <- years[tidx]
  mons <- months[tidx]
  nt   <- dim(arr)[3]
  nlon <- length(lons)
  nlat <- length(lats)

  # Monthly climatology over clim.range
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

  # [time x cells]
  mat  <- matrix(anom, nrow = nt)

  # Keep cells with data present in all time steps.
  # Relax to 95% if needed (e.g. sparse sea-ice months).
  complete <- colSums(!is.na(mat))
  good <- which(complete == nt)
  if (length(good) < 50) good <- which(complete >= nt * 0.95)

  mat2 <- mat[, good]

  # Fill any residual sparse NAs with the cell temporal mean (= ~0 for anomalies)
  for (j in seq_len(ncol(mat2))) {
    bad <- is.na(mat2[, j])
    if (any(bad)) mat2[bad, j] <- mean(mat2[, j], na.rm = TRUE)
  }

  # Detrend each cell
  t.seq <- seq_len(nt)
  mat2  <- apply(mat2, 2, function(x) residuals(lm(x ~ t.seq)))

  # Latitude weights: sqrt(cos(lat))
  grid    <- expand.grid(lon = lons, lat = lats)
  weights <- sqrt(cos(grid$lat[good] * pi / 180))

  # Weighted data matrix [time x cells]
  mat.w <- sweep(mat2, 2, weights, "*")

  # Temporal covariance [time x time] — small matrix, numerically stable
  C   <- tcrossprod(mat.w) / (nt - 1)
  eig <- eigen(C, symmetric = TRUE)

  # Leading PC (time scores)
  pc1 <- eig$vectors[, 1]

  # Recover EOF1 spatial loadings by projecting weighted matrix onto PC1
  eof1.good <- as.vector(crossprod(mat.w, pc1))
  eof1.good <- eof1.good / sqrt(sum(eof1.good^2))

  loading       <- rep(NA_real_, nlon * nlat)
  loading[good] <- eof1.good

  df         <- expand.grid(lon = lons, lat = lats)
  df$loading <- loading
  df
}

# Sign convention: flip so dominant lobe is negative (PDO convention: cool central NP)
orient_eof <- function(df) {
  if (mean(df$loading, na.rm = TRUE) > 0) df$loading <- -df$loading
  df
}

# Coarsen a loading data frame to ~res degree bins for faster plotting
coarsen_df <- function(df, res = 1) {
  df %>%
    mutate(lon = round(lon / res) * res,
           lat = round(lat / res) * res) %>%
    group_by(lon, lat) %>%
    summarise(loading = mean(loading, na.rm = TRUE), .groups = "drop")
}

# ---- PLOTTING ----

eof_panel <- function(df, title, col.lim) {
  df <- coarsen_df(df, res = 2)
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

# ---- LOAD CESM2 MEMBER ----
# Crops to the common domain (lon.min360-lon.max360, lat.min-lat.max) on load.
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

# Enforce common domain (SST file already covers 110-250E but confirm)
lon.idx.sst <- which(to360(lons.sst) >= lon.min360 & to360(lons.sst) <= lon.max360)
lat.idx.sst <- which(lats.sst >= lat.min & lats.sst <= lat.max)
sst.crop    <- sst.raw[lon.idx.sst, lat.idx.sst, ]
lons.sst2   <- lons.sst[lon.idx.sst]
lats.sst2   <- lats.sst[lat.idx.sst]
rm(sst.raw); gc()

message("Computing ERA5 SST EOF1...")
eof.sst.obs <- orient_eof(
  compute_eof1(sst.crop, lons.sst2, lats.sst2, yrs.e, mons.e)
)
rm(sst.crop); gc()

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

# Crop to same domain as SST and CESM2: 110-250E, 20-66N
lon.idx.slp <- which(to360(lons.slp) >= lon.min360 & to360(lons.slp) <= lon.max360)
lat.idx.slp <- which(lats.slp >= lat.min & lats.slp <= lat.max)
slp.crop    <- slp.raw[lon.idx.slp, lat.idx.slp, ]
lons.slp2   <- lons.slp[lon.idx.slp]
lats.slp2   <- lats.slp[lat.idx.slp]
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

save_panel <- function(p, filename) {
  ggsave(filename, plot = p, width = 16, height = 10, dpi = 150)
  message("Saved: ", filename)
}

save_panel(make_12panel(eof.sst.obs, fcm.sst.eofs, "SST", "FCM"),
           "./Figures/EOF1_SST_FCM.png")
save_panel(make_12panel(eof.sst.obs, mdm.sst.eofs, "SST", "MDM"),
           "./Figures/EOF1_SST_MDM.png")
save_panel(make_12panel(eof.slp.obs, fcm.slp.eofs, "SLP", "FCM"),
           "./Figures/EOF1_SLP_FCM.png")
save_panel(make_12panel(eof.slp.obs, mdm.slp.eofs, "SLP", "MDM"),
           "./Figures/EOF1_SLP_MDM.png")
