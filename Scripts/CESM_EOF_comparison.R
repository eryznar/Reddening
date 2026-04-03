# PURPOSE: Compare ERA5 and CESM2 FCM/MDM EOF1 loadings for SST and SLP
# Author: Mike Litzow
#
# Produces four 4x3 panel plots (12 panels each):
#   Panel 1    = ERA5 observations
#   Panels 2-12 = 11 CESM2 ensemble members
#
# Plot 1: SST EOF1 — ERA5 obs vs FCM members
# Plot 2: SST EOF1 — ERA5 obs vs MDM members
# Plot 3: SLP EOF1 — ERA5 obs vs FCM members
# Plot 4: SLP EOF1 — ERA5 obs vs MDM members
#
# Common period: 1950-2014 (CESM2 ends 2014)
# Season: winter (Nov-Mar); year assigned to January year
# Domain: 20-66N, 110-250E

source("./Scripts/load.libs.functions.R")

# ---- SETTINGS ----
win.months <- c(11, 12, 1, 2, 3)
yr.range   <- c(1950, 2014)
n.members  <- 11
lat.min    <- 20;  lat.max    <- 66
lon.min360 <- 110; lon.max360 <- 250

mapWorld <- map_data("world", wrap = c(20, 380))

# ---- HELPER FUNCTIONS ----

# Compute winter-mean anomalies from monthly [lon x lat x time] array.
# Climatology computed over all months in the input array.
# Returns list: arr [lon x lat x n.winter.years], years vector.
winter_anom <- function(arr, years, months) {
  clim <- array(NA, dim = c(dim(arr)[1], dim(arr)[2], 12))
  for (m in 1:12) {
    mi <- which(months == m)
    if (length(mi) > 0)
      clim[, , m] <- apply(arr[, , mi, drop = FALSE], c(1, 2),
                           mean, na.rm = TRUE)
  }
  anom <- arr
  for (t in seq_along(months))
    anom[, , t] <- arr[, , t] - clim[, , months[t]]

  wi   <- which(months %in% win.months)
  wyrs <- ifelse(months[wi] %in% c(11, 12), years[wi] + 1L, years[wi])
  uyrs <- sort(unique(wyrs[wyrs >= yr.range[1] & wyrs <= yr.range[2]]))

  out <- array(NA, dim = c(dim(arr)[1], dim(arr)[2], length(uyrs)))
  for (i in seq_along(uyrs)) {
    ti       <- wi[wyrs == uyrs[i]]
    out[,,i] <- apply(anom[, , ti, drop = FALSE], c(1, 2), mean, na.rm = TRUE)
  }
  list(arr = out, years = uyrs)
}

# Compute EOF1 loadings from [lon x lat x time] array.
# Returns data frame: lon, lat, loading.
compute_eof1 <- function(arr, lons, lats) {
  nlon <- length(lons); nlat <- length(lats); nt <- dim(arr)[3]
  mat  <- matrix(arr, nrow = nt)
  good <- which(colSums(!is.na(mat)) > nt * 0.5)
  mat2 <- mat[, good]
  mat2 <- apply(mat2, 2, function(x) {
    ok <- !is.na(x)
    if (sum(ok) > 5) x[ok] <- residuals(lm(x[ok] ~ seq_len(sum(ok))))
    x
  })
  grid    <- expand.grid(lon = lons, lat = lats)
  weights <- sqrt(cos(grid$lat[good] * pi / 180))
  pca     <- svd.triplet(cov(mat2, use = "pairwise.complete.obs"),
                         col.w = weights)
  loadings         <- rep(NA_real_, nlon * nlat)
  loadings[good]   <- pca$U[, 1]
  df               <- expand.grid(lon = lons, lat = lats)
  df$loading       <- loadings
  df
}

# Sign convention: flip so overall mean loading is negative
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
      plot.title      = element_text(size = 7, hjust = 0.5),
      axis.title      = element_blank(),
      axis.text       = element_text(size = 5),
      panel.grid      = element_blank(),
      legend.position = "right",
      legend.key.height = unit(0.8, "cm")
    )
}

# 4x3 panel plot: obs + n.members model members
make_12panel <- function(obs.df, member.dfs, var.label, model.label) {
  all.load <- c(obs.df$loading, unlist(lapply(member.dfs, `[[`, "loading")))
  col.lim  <- max(abs(all.load), na.rm = TRUE)

  panels <- c(
    list(eof_panel(obs.df, paste0("ERA5 Obs (", var.label, ")"), col.lim)),
    lapply(seq_along(member.dfs), function(i)
      eof_panel(member.dfs[[i]],
                paste0(model.label, " member ", i), col.lim))
  )

  wrap_plots(panels, ncol = 4) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste(
        "EOF1 Loadings — ERA5 vs CESM2", model.label, "(", var.label, ")",
        "\nWinter (Nov-Mar), 1950-2014"
      ),
      theme = theme(plot.title = element_text(size = 10, hjust = 0.5))
    )
}

# ---- LOAD CESM2 MEMBER ----
# Returns list: arr [lon x lat x time], lons, lats, years, months
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

  nt     <- dim(dat)[3]
  dates  <- seq(as.Date("1850-01-15"), by = "month", length.out = nt)
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

dates.e  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
yrs.e    <- as.integer(format(dates.e, "%Y"))
mons.e   <- as.integer(format(dates.e, "%m"))

idx.e    <- which(yrs.e >= yr.range[1] & yrs.e <= yr.range[2])
wa.sst   <- winter_anom(sst.raw[, , idx.e], yrs.e[idx.e], mons.e[idx.e])
eof.sst.obs <- orient_eof(compute_eof1(wa.sst$arr, lons.sst, lats.sst))

# ---- ERA5 SLP ----
message("Loading ERA5 SLP...")
nc       <- nc_open("./Data/era5_slp_NP_1950_2025.nc")
slp.raw  <- ncvar_get(nc, "msl") / 100     # Pa -> hPa
lons.slp <- ncvar_get(nc, "longitude")
lats.slp <- ncvar_get(nc, "latitude")
time     <- ncvar_get(nc, "valid_time")
nc_close(nc)

dates.s  <- as.Date(as.POSIXct(time, origin = "1970-01-01", tz = "UTC"))
yrs.s    <- as.integer(format(dates.s, "%Y"))
mons.s   <- as.integer(format(dates.s, "%m"))

# Restrict to common lat/lon domain (SLP file has broader coverage)
lat.idx.slp <- which(lats.slp >= lat.min & lats.slp <= lat.max)
lon.idx.slp <- which(to360(lons.slp) >= lon.min360 &
                       to360(lons.slp) <= lon.max360)
slp.crop  <- slp.raw[lon.idx.slp, lat.idx.slp, ]
lons.slp2 <- lons.slp[lon.idx.slp]
lats.slp2 <- lats.slp[lat.idx.slp]

idx.s    <- which(yrs.s >= yr.range[1] & yrs.s <= yr.range[2])
wa.slp   <- winter_anom(slp.crop[, , idx.s], yrs.s[idx.s], mons.s[idx.s])
eof.slp.obs <- orient_eof(compute_eof1(wa.slp$arr, lons.slp2, lats.slp2))

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

# ---- PROCESS CESM2 MEMBERS ----

eof_cesm <- function(file, varname, scale.pa = FALSE) {
  d <- load_cesm(file, varname)
  if (scale.pa) d$arr <- d$arr / 100
  wa <- winter_anom(d$arr, d$years, d$months)
  orient_eof(compute_eof1(wa$arr, d$lons, d$lats))
}

message("Processing FCM SST members (", length(fcm.sst.files), ")...")
fcm.sst.eofs <- lapply(seq_along(fcm.sst.files), function(i) {
  message("  member ", i)
  eof_cesm(fcm.sst.files[i], "SST")
})

message("Processing MDM SST members (", length(mdm.sst.files), ")...")
mdm.sst.eofs <- lapply(seq_along(mdm.sst.files), function(i) {
  message("  member ", i)
  eof_cesm(mdm.sst.files[i], "SST")
})

message("Processing FCM SLP members (", length(fcm.slp.files), ")...")
fcm.slp.eofs <- lapply(seq_along(fcm.slp.files), function(i) {
  message("  member ", i)
  eof_cesm(fcm.slp.files[i], "PSL", scale.pa = TRUE)
})

message("Processing MDM SLP members (", length(mdm.slp.files), ")...")
mdm.slp.eofs <- lapply(seq_along(mdm.slp.files), function(i) {
  message("  member ", i)
  eof_cesm(mdm.slp.files[i], "PSL", scale.pa = TRUE)
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
