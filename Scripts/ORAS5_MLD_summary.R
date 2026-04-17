# PURPOSE: Summarize ORAS5 MLD output — cell-wise mean maps
# Author: Mike Litzow
#
# Produces two maps saved to Figures/:
#   1. Annual mean MLD  (all months, all years)
#   2. Winter mean MLD  (November-March, all years)

source("./Scripts/load.libs.functions.R")

# ---- PRECOMPUTE MEANS VIA CDO (avoids loading full 2GB array into R) ----
# Data is on the ORCA025 tripolar grid (~1/4 degree resolution).
# CDO streams the file without loading all time steps at once.

cdo <- "/usr/local/bin/cdo"
src <- "./Data/oras5_mld_NP_1958_2025.nc"

if (!file.exists("./Data/oras5_mld_NP_mean_all.nc")) {
  message("Computing all-month mean via CDO...")
  system(paste(cdo, "timmean", src, "./Data/oras5_mld_NP_mean_all.nc"))
}

if (!file.exists("./Data/oras5_mld_NP_mean_win.nc")) {
  message("Computing winter (Nov-Mar) mean via CDO...")
  system(paste(cdo, "timmean -selmon,11,12,1,2,3", src,
               "./Data/oras5_mld_NP_mean_win.nc"))
}

# ---- LOAD PRECOMPUTED MEANS (small files, single time step each) ----
load_mean_nc <- function(file) {
  nc   <- nc_open(file)
  raw  <- ncvar_get(nc, "somxl030")
  mld  <- if (length(dim(raw)) == 3) raw[, , 1] else raw
  lons <- ncvar_get(nc, "nav_lon")
  lats <- ncvar_get(nc, "nav_lat")
  nc_close(nc)
  mld[mld > 1e10] <- NA
  list(mld = mld, lons = lons, lats = lats)
}

d_all <- load_mean_nc("./Data/oras5_mld_NP_mean_all.nc")
d_win <- load_mean_nc("./Data/oras5_mld_NP_mean_win.nc")

mean_all <- d_all$mld;  lons <- d_all$lons;  lats <- d_all$lats
mean_win <- d_win$mld

# ---- BUILD DATA FRAMES ----
# Flatten 2D grid
df_all <- data.frame(
  lon = as.vector(lons),
  lat = as.vector(lats),
  mld = as.vector(mean_all)
) %>%
  filter(!is.na(mld), mld > 0) %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon)) %>%
  filter(lat >= 20, lat <= 66, lon >= 110, lon <= 250)

df_win <- data.frame(
  lon = as.vector(lons),
  lat = as.vector(lats),
  mld = as.vector(mean_win)
) %>%
  filter(!is.na(mld), mld > 0) %>%
  mutate(lon = ifelse(lon < 0, lon + 360, lon)) %>%
  filter(lat >= 20, lat <= 66, lon >= 110, lon <= 250)

# ---- MAP SETUP ----
mapWorld <- map_data("world", wrap = c(20, 380))

# Study region polygons in 0-360 longitude space
goa.poly <- data.frame(
  lon = c(198, 198, 203, 205, 208, 225, 231, 201),
  lat = c(54, 55.5, 57.5, 59, 61, 61, 54, 54),
  region = "GOA"
)
ebs.poly <- data.frame(
  lon = c(183, 183, 203, 203, 191),
  lat = c(53, 65, 65, 57.5, 53),
  region = "EBS"
)
polys <- bind_rows(goa.poly, ebs.poly)

mld_map <- function(df, title) {
  ggplot() +
    geom_point(data = df, aes(x = lon, y = lat, color = mld),
               shape = 15, size = 0.4) +
    geom_polygon(data = mapWorld,
                 aes(x = long, y = lat, group = group),
                 fill = "gray80", color = "gray50", linewidth = 0.2) +
    geom_polygon(data = polys,
                 aes(x = lon, y = lat, group = region),
                 fill = NA, color = "black", linewidth = 0.7) +
    scale_color_gradientn(
      colors = oceColorsPalette(64),
      name   = "MLD (m)"
    ) +
    coord_cartesian(xlim = c(110, 250), ylim = c(20, 66)) +
    labs(title = title,
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(
      plot.title      = element_text(hjust = 0.5),
      legend.position = "right"
    )
}

# ---- PLOTS ----
dir.create("./Figures", showWarnings = FALSE)

p_all <- mld_map(df_all,
  "ORAS5 MLD — Cell-wise Mean, All Months (1958-2025)")

p_win <- mld_map(df_win,
  "ORAS5 MLD — Cell-wise Mean, Winter (Nov-Mar, 1958-2025)")

ggsave("./Figures/ORAS5_MLD_mean_all_months.png",
       plot = p_all, width = 10, height = 5, dpi = 150)
ggsave("./Figures/ORAS5_MLD_mean_winter.png",
       plot = p_win, width = 10, height = 5, dpi = 150)

message("Saved: Figures/ORAS5_MLD_mean_all_months.png")
message("Saved: Figures/ORAS5_MLD_mean_winter.png")

# ---- TIME SERIES — GOA AND EBS ----
# Strategy: use CDO sellonlatbox to extract small bounding boxes (manageable in R),
# then apply polygon masking inside R before computing weighted spatial means.
#
# Polygons are defined in 0-360 space; ORAS5 nav_lon is -180/180.
# Bounding boxes converted accordingly:
#   GOA: lon 198-231E = -162 to -129W, lat 54-61N
#   EBS: lon 183-203E = -177 to -157W, lat 53-65N

# CDO bounding box extraction (run once, cached)
goa.box.file <- "./Data/oras5_mld_GOA_box.nc"
ebs.box.file <- "./Data/oras5_mld_EBS_box.nc"

if (!file.exists(goa.box.file)) {
  message("Extracting GOA bounding box via CDO...")
  system(paste(cdo, "sellonlatbox,-162,-129,54,61", src, goa.box.file))
}
if (!file.exists(ebs.box.file)) {
  message("Extracting EBS bounding box via CDO...")
  system(paste(cdo, "sellonlatbox,-177,-157,53,65", src, ebs.box.file))
}

# Load a bounding-box file, return list with mld [x,y,t], nav_lon, nav_lat, dates
load_box <- function(file) {
  nc     <- nc_open(file)
  mld    <- ncvar_get(nc, "somxl030")
  lons   <- ncvar_get(nc, "nav_lon")
  lats   <- ncvar_get(nc, "nav_lat")
  time   <- ncvar_get(nc, "time_counter")
  tunits <- ncatt_get(nc, "time_counter", "units")$value
  nc_close(nc)
  mld[mld > 1e10] <- NA
  origin <- sub(".*since ", "", tunits)
  dates  <- as.Date(as.POSIXct(time, origin = origin, tz = "UTC"))
  list(mld = mld, lons = lons, lats = lats, dates = dates)
}

goa.d <- load_box(goa.box.file)
ebs.d <- load_box(ebs.box.file)

# Build sf polygon objects for masking (coordinates in -180/180)
make_sf_poly <- function(lon360, lat) {
  lon180 <- ifelse(lon360 > 180, lon360 - 360, lon360)
  coords  <- rbind(cbind(lon180, lat), c(lon180[1], lat[1]))
  st_sf(geometry = st_sfc(st_polygon(list(coords)), crs = 4326))
}

goa.sf <- make_sf_poly(
  c(198, 198, 203, 205, 208, 225, 231, 201),
  c(54,  55.5, 57.5, 59,  61,  61,  54,  54)
)
ebs.sf <- make_sf_poly(
  c(183, 183, 203, 203, 191),
  c(53,  65,  65,  57.5, 53)
)

# Compute latitude-weighted spatial mean time series for one region
region_ts <- function(d, poly.sf) {
  lons <- as.vector(d$lons)
  lats <- as.vector(d$lats)

  # Point-in-polygon mask
  grid.sf <- st_as_sf(data.frame(lon = lons, lat = lats),
                      coords = c("lon", "lat"), crs = 4326)
  in.poly <- lengths(st_within(grid.sf, poly.sf)) > 0

  if (sum(in.poly) == 0) stop("No grid cells inside polygon")

  weights <- cos(lats[in.poly] * pi / 180)
  nt      <- dim(d$mld)[3]

  # Flatten spatial dims: mld is [x, y, t] or [xy, t] after matrix()
  mat <- matrix(d$mld, nrow = prod(dim(d$mld)[1:2]), ncol = nt)

  sapply(seq_len(nt), function(t) {
    vals <- mat[in.poly, t]
    if (all(is.na(vals))) return(NA_real_)
    weighted.mean(vals, w = weights, na.rm = TRUE)
  })
}

message("Computing GOA spatial mean time series...")
goa.ts <- region_ts(goa.d, goa.sf)
message("Computing EBS spatial mean time series...")
ebs.ts <- region_ts(ebs.d, ebs.sf)

dates  <- goa.d$dates
years  <- as.integer(format(dates, "%Y"))
months <- as.integer(format(dates, "%m"))

# Build monthly data frame, compute anomaly and detrended anomaly
make_monthly <- function(ts.vec, yrs, mons) {
  df <- data.frame(year = yrs, month = mons, mld = ts.vec)
  clim <- df %>%
    group_by(month) %>%
    summarise(clim = mean(mld, na.rm = TRUE), .groups = "drop")
  df <- df %>%
    left_join(clim, by = "month") %>%
    mutate(anom      = mld - clim,
           anom.detr = residuals(lm(anom ~ seq_along(anom))))
  df
}

goa.monthly <- make_monthly(goa.ts, years, months)
ebs.monthly <- make_monthly(ebs.ts, years, months)

# Winter means (Nov-Mar, year = January year)
make_winter <- function(df, region) {
  df %>%
    filter(month %in% c(11, 12, 1, 2, 3)) %>%
    mutate(win.year = ifelse(month %in% c(11, 12), year + 1L, year)) %>%
    group_by(win.year) %>%
    summarise(raw       = mean(mld,       na.rm = TRUE),
              anom      = mean(anom,      na.rm = TRUE),
              anom.detr = mean(anom.detr, na.rm = TRUE),
              .groups   = "drop") %>%
    rename(year = win.year) %>%
    mutate(region = region)
}

winter <- bind_rows(
  make_winter(goa.monthly, "GOA"),
  make_winter(ebs.monthly, "EBS")
)

region.colors <- c("GOA" = "steelblue4", "EBS" = "darkorange4")

ts_long <- winter %>%
  pivot_longer(c(raw, anom, anom.detr),
               names_to = "series", values_to = "value") %>%
  mutate(series = factor(series,
    levels = c("raw", "anom", "anom.detr"),
    labels = c("Raw MLD (m)",
               "Monthly anomaly (m)",
               "Detrended anomaly (m)")))

p_ts <- ggplot(ts_long, aes(x = year, y = value, color = region)) +
  geom_line() +
  geom_point(size = 1.2) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "gray50") +
  scale_color_manual(values = region.colors) +
  facet_wrap(~ series, ncol = 1, scales = "free_y") +
  labs(title = "ORAS5 MLD — GOA and EBS Winter (Nov-Mar) Time Series",
       x = "Year", y = NULL, color = "Region") +
  theme_bw() +
  theme(plot.title       = element_text(hjust = 0.5),
        strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

ggsave("./Figures/ORAS5_MLD_timeseries_winter.png",
       plot = p_ts, width = 10, height = 8, dpi = 150)
message("Saved: Figures/ORAS5_MLD_timeseries_winter.png")

# ---- CELLWISE REGRESSION: SLP PC1 -> MLD ----
# Regress winter-mean detrended MLD anomaly (y) on winter SLP PC1 (x)
# at each grid cell. Restricted to NE Pacific: 170-250E, >=30N.
# Uses GLS with AR(1) residuals (nlme::gls + corAR1).
# Contour line drawn around cells with p < 0.05.
#
# CDO pipeline (run once, cached):
#   1. ymonmean   -> monthly climatology
#   2. ymonsub    -> monthly anomalies
#   3. detrend    -> detrended monthly anomalies
#   4. selmon     -> keep Nov-Mar only

mld.winmon.file <- "./Data/oras5_mld_NP_detr_anom_winmon.nc"

if (!file.exists(mld.winmon.file)) {
  message("CDO: computing detrended anomalies and selecting winter months...")
  system(paste(
    cdo, "selmon,11,12,1,2,3",
    "-detrend",
    "-ymonsub", src,
    paste0("-ymonmean ", src),
    mld.winmon.file
  ))
}

# GLS with AR(1) residuals for each cell — skip everything if cached results exist
reg.cache <- "./Output/oras5_mld_slp_pc1_regression.csv"
if (file.exists(reg.cache)) {
  message("Loading cached regression results...")
  reg.df <- read.csv(reg.cache)
} else {
  # Pre-crop to NE Pacific bounding box via CDO before loading into R
  mld.ne.file <- "./Data/oras5_mld_NEP_detr_anom_winmon.nc"
  if (!file.exists(mld.ne.file)) {
    message("CDO: cropping to NE Pacific bounding box...")
    system(paste(cdo, "sellonlatbox,-190,-110,30,66",
                 mld.winmon.file, mld.ne.file))
  }

  message("Loading cellwise winter-month MLD anomalies (NE Pacific)...")
  nc     <- nc_open(mld.ne.file)
  mld.w  <- ncvar_get(nc, "somxl030")
  lons.w <- ncvar_get(nc, "nav_lon")
  lats.w <- ncvar_get(nc, "nav_lat")
  time.w <- ncvar_get(nc, "time_counter")
  tun.w  <- ncatt_get(nc, "time_counter", "units")$value
  nc_close(nc)

  mld.w[mld.w > 1e10] <- NA

  origin.w   <- sub(".*since ", "", tun.w)
  dates.w    <- as.Date(as.POSIXct(time.w, origin = origin.w, tz = "UTC"))
  years.w    <- as.integer(format(dates.w, "%Y"))
  months.w   <- as.integer(format(dates.w, "%m"))
  winyears.w <- ifelse(months.w %in% c(11, 12), years.w + 1L, years.w)

  lons.w360  <- ifelse(lons.w < 0, lons.w + 360, lons.w)
  ne.mask    <- lons.w360 >= 170 & lons.w360 <= 250 & lats.w >= 30

  pc1        <- read.csv("./Output/NP_PC1_winter.csv")
  common.yrs <- sort(intersect(unique(winyears.w), pc1$year))
  pc1.vec    <- pc1$PC1[match(common.yrs, pc1$year)]
  n.yrs      <- length(common.yrs)

  nx   <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
  mat  <- matrix(mld.w, nrow = nx * ny, ncol = nt)
  good.cells <- which(ne.mask & rowSums(!is.na(mat)) > nt * 0.5)
  rm(mld.w); gc()

  message("Computing cellwise winter means...")
  win.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = n.yrs)
  for (k in seq_along(common.yrs)) {
    ti <- which(winyears.w == common.yrs[k])
    win.mat[, k] <- rowMeans(mat[good.cells, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat); gc()

  message("Fitting GLS AR(1) regressions (this may take several minutes)...")
beta.vec <- rep(NA_real_, length(good.cells))
pval.vec <- rep(NA_real_, length(good.cells))

df.reg <- data.frame(y = NA_real_, x = pc1.vec)

for (k in seq_along(good.cells)) {
  y <- win.mat[k, ]
  if (sum(!is.na(y)) < 10) next
  df.reg$y <- y
  fit <- tryCatch(
    gls(y ~ x, data = df.reg,
        correlation = corAR1(form = ~ 1),
        method = "ML", na.action = na.omit),
    error = function(e) NULL
  )
  if (is.null(fit)) next
  s <- summary(fit)$tTable
  if (nrow(s) < 2) next
  beta.vec[k] <- s[2, 1]
  pval.vec[k] <- s[2, 4]
}

# Assemble full result grid
beta.full <- rep(NA_real_, nx * ny)
pval.full <- rep(NA_real_, nx * ny)
beta.full[good.cells] <- beta.vec
pval.full[good.cells] <- pval.vec

reg.df <- data.frame(
  lon  = as.vector(lons.w360),
  lat  = as.vector(lats.w),
  beta = beta.full,
  pval = pval.full
) %>%
  filter(!is.na(beta)) %>%
  filter(lon >= 170, lon <= 250, lat >= 30)

# Cache regression results to avoid recomputing
  write.csv(reg.df, "./Output/oras5_mld_slp_pc1_regression.csv", row.names = FALSE)
  message("Regression results saved to Output/oras5_mld_slp_pc1_regression.csv")
} # end else (regression computed fresh)

# common.yrs is only defined above when the cache is absent (fresh computation).
# When loading from cache, derive it from the PC1 file so the subtitle renders.
if (!exists("common.yrs")) {
  pc1.tmp    <- read.csv("./Output/NP_PC1_winter.csv")
  common.yrs <- pc1.tmp$year
  rm(pc1.tmp)
}

col.lim <- max(abs(reg.df$beta), na.rm = TRUE)

# Significance contour: bin to regular grid then geom_contour
reg.grid <- reg.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            pval = mean(pval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(170, 250, by = 0.5),
           lat = seq(30,  66,  by = 0.5))

p_reg <- ggplot() +
  geom_point(data = filter(reg.df, !is.na(beta)),
             aes(x = lon, y = lat, color = beta),
             shape = 15, size = 0.5) +
  geom_contour(data = filter(reg.grid, !is.na(pval)),
               aes(x = lon, y = lat, z = pval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region),
               fill = NA, color = "black", linewidth = 0.7) +
  scale_color_gradientn(
    colors   = rev(oceColorsPalette(64)),
    limits   = c(-col.lim, col.lim),
    values   = scales::rescale(c(-col.lim, 0, col.lim)),
    na.value = "white",
    name     = "β (m / PC1)"
  ) +
  coord_cartesian(xlim = c(170, 250), ylim = c(30, 66)) +
  labs(
    title    = "Regression: SLP PC1 → Cellwise Winter MLD Anomaly (detrended)",
    subtitle = paste0("GLS AR(1); contour = p < 0.05; ",
                      min(common.yrs), "–", max(common.yrs)),
    x = "Longitude", y = "Latitude"
  ) +
  theme_bw() +
  theme(plot.title    = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        legend.position = "right")

ggsave("./Figures/ORAS5_MLD_SLP_PC1_regression.png",
       plot = p_reg, width = 9, height = 5, dpi = 150)
message("Saved: Figures/ORAS5_MLD_SLP_PC1_regression.png")

# ---- ERA-SPECIFIC REGRESSIONS: AL low-variability vs high-variability ----
# Low-variability AL era:  1998-2010
# High-variability AL era: 2016-2025
# Uses same win.mat (cellwise winter means) and pc1 as the full-period regression.
# win.mat columns correspond to common.yrs — subset by year range.

era_regression <- function(yr.min, yr.max, cache.file) {
  if (file.exists(cache.file)) {
    message("Loading cached era regression: ", cache.file)
    return(read.csv(cache.file))
  }

  era.idx <- which(common.yrs >= yr.min & common.yrs <= yr.max)
  if (length(era.idx) < 5) stop("Too few years in era: ", yr.min, "-", yr.max)

  y.mat   <- win.mat[, era.idx]
  x.vec   <- pc1.vec[era.idx]
  n.cells <- nrow(y.mat)

  beta.v <- rep(NA_real_, n.cells)
  pval.v <- rep(NA_real_, n.cells)
  df.r   <- data.frame(y = NA_real_, x = x.vec)

  message("Fitting GLS AR(1) regressions for ", yr.min, "-", yr.max, "...")
  for (k in seq_len(n.cells)) {
    y <- y.mat[k, ]
    if (sum(!is.na(y)) < 5) next
    df.r$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.r,
          correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.v[k] <- s[2, 1]
    pval.v[k] <- s[2, 4]
  }

  beta.full <- rep(NA_real_, nx * ny)
  pval.full <- rep(NA_real_, nx * ny)
  beta.full[good.cells] <- beta.v
  pval.full[good.cells] <- pval.v

  out <- data.frame(
    lon  = as.vector(lons.w360),
    lat  = as.vector(lats.w),
    beta = beta.full,
    pval = pval.full
  ) %>%
    filter(!is.na(beta)) %>%
    filter(lon >= 170, lon <= 250, lat >= 30)

  write.csv(out, cache.file, row.names = FALSE)
  message("Saved: ", cache.file)
  out
}

# Run or load era regressions
# Note: win.mat, good.cells, lons.w360, lats.w, nx, ny, common.yrs, pc1.vec
# must be in environment. If only the full-period cache was loaded above,
# we need to reload the NE Pacific array to recompute win.mat.
era.cache.low  <- "./Output/oras5_mld_slp_pc1_regression_low_era.csv"
era.cache.high <- "./Output/oras5_mld_slp_pc1_regression_high_era.csv"

if (!file.exists(era.cache.low) || !file.exists(era.cache.high)) {
  if (!exists("win.mat")) {
    message("Reloading NE Pacific MLD array for era regressions...")
    mld.ne.file <- "./Data/oras5_mld_NEP_detr_anom_winmon.nc"
    nc     <- nc_open(mld.ne.file)
    mld.w  <- ncvar_get(nc, "somxl030")
    lons.w <- ncvar_get(nc, "nav_lon")
    lats.w <- ncvar_get(nc, "nav_lat")
    time.w <- ncvar_get(nc, "time_counter")
    tun.w  <- ncatt_get(nc, "time_counter", "units")$value
    nc_close(nc)
    mld.w[mld.w > 1e10] <- NA
    origin.w   <- sub(".*since ", "", tun.w)
    dates.w    <- as.Date(as.POSIXct(time.w, origin = origin.w, tz = "UTC"))
    years.w    <- as.integer(format(dates.w, "%Y"))
    months.w   <- as.integer(format(dates.w, "%m"))
    winyears.w <- ifelse(months.w %in% c(11, 12), years.w + 1L, years.w)
    lons.w360  <- ifelse(lons.w < 0, lons.w + 360, lons.w)
    ne.mask    <- lons.w360 >= 170 & lons.w360 <= 250 & lats.w >= 30
    pc1        <- read.csv("./Output/NP_PC1_winter.csv")
    common.yrs <- sort(intersect(unique(winyears.w), pc1$year))
    pc1.vec    <- pc1$PC1[match(common.yrs, pc1$year)]
    n.yrs      <- length(common.yrs)
    nx <- dim(mld.w)[1]; ny <- dim(mld.w)[2]; nt <- dim(mld.w)[3]
    mat <- matrix(mld.w, nrow = nx * ny, ncol = nt)
    good.cells <- which(ne.mask & rowSums(!is.na(mat)) > nt * 0.5)
    rm(mld.w); gc()
    win.mat <- matrix(NA_real_, nrow = length(good.cells), ncol = n.yrs)
    for (k in seq_along(common.yrs)) {
      ti <- which(winyears.w == common.yrs[k])
      win.mat[, k] <- rowMeans(mat[good.cells, ti, drop = FALSE], na.rm = TRUE)
    }
    rm(mat); gc()
  }
}

reg.low  <- era_regression(1998, 2010, era.cache.low)
reg.high <- era_regression(2016, 2025, era.cache.high)

# Shared color limit across both panels
col.lim.era <- max(abs(c(reg.low$beta, reg.high$beta)), na.rm = TRUE)

make_reg_panel <- function(df, title.label) {
  grid.df <- df %>%
    mutate(lon = round(lon / 0.5) * 0.5,
           lat = round(lat / 0.5) * 0.5) %>%
    group_by(lon, lat) %>%
    summarise(beta = mean(beta, na.rm = TRUE),
              pval = mean(pval, na.rm = TRUE),
              .groups = "drop") %>%
    complete(lon = seq(170, 250, by = 0.5),
             lat = seq(30,  66,  by = 0.5))

  ggplot() +
    geom_point(data = filter(df, !is.na(beta)),
               aes(x = lon, y = lat, color = beta),
               shape = 15, size = 0.5) +
    geom_contour(data = filter(grid.df, !is.na(pval)),
                 aes(x = lon, y = lat, z = pval),
                 breaks = 0.05, color = "black", linewidth = 0.4) +
    geom_polygon(data = mapWorld,
                 aes(x = long, y = lat, group = group),
                 fill = "gray80", color = "gray50", linewidth = 0.2) +
    geom_polygon(data = polys,
                 aes(x = lon, y = lat, group = region),
                 fill = NA, color = "black", linewidth = 0.7) +
    scale_color_gradientn(
      colors   = rev(oceColorsPalette(64)),
      limits   = c(-col.lim.era, col.lim.era),
      values   = scales::rescale(c(-col.lim.era, 0, col.lim.era)),
      na.value = "white",
      name     = "β (m / PC1)"
    ) +
    coord_cartesian(xlim = c(170, 250), ylim = c(30, 66)) +
    labs(title = title.label,
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(plot.title      = element_text(hjust = 0.5),
          legend.position = "right")
}

p.low  <- make_reg_panel(reg.low,  "Low AL variability era (1998-2010)")
p.high <- make_reg_panel(reg.high, "High AL variability era (2016-2025)")

p.era <- p.low / p.high +
  plot_annotation(
    title    = "Regression by Era: SLP PC1 → Cellwise Winter MLD Anomaly (detrended)",
    subtitle = "GLS AR(1); contour = p < 0.05",
    theme    = theme(plot.title    = element_text(hjust = 0.5),
                     plot.subtitle = element_text(hjust = 0.5, size = 9))
  )

ggsave("./Figures/ORAS5_MLD_SLP_PC1_regression_by_era.png",
       plot = p.era, width = 9, height = 8, dpi = 150)
message("Saved: Figures/ORAS5_MLD_SLP_PC1_regression_by_era.png")

# ---- 15-YR ROLLING AL SD vs MLD AR(1) — dual y-axis, GOA and EBS panels ----
# Uses detrended winter MLD anomaly for AR(1); AL SLP SD from pre-computed CSV.

roll_ar1_mld <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right",
            FUN = function(v) acf(v, lag.max = 1, plot = FALSE)$acf[2])
}
roll_sd_al <- function(x, width = 15) {
  rollapply(x, width = width, fill = NA, align = "right", FUN = sd)
}

# Build wide winter MLD frame from the monthly regional data frames
goa.win <- make_winter(goa.monthly, "GOA") %>% select(year, GOA = anom.detr)
ebs.win <- make_winter(ebs.monthly, "EBS") %>% select(year, EBS = anom.detr)
winter.wide <- inner_join(goa.win, ebs.win, by = "year") %>%
  arrange(year) %>%
  mutate(GOA.ar1 = roll_ar1_mld(GOA),
         EBS.ar1 = roll_ar1_mld(EBS))

# Load AL SLP anomaly; compute 15-yr rolling SD, z-scored
al.slp <- read.csv("./Output/AL_winter_SLP_anomaly.csv") %>%
  arrange(year) %>%
  mutate(AL.sd = roll_sd_al(SLP),
         AL.sd = as.numeric(scale(AL.sd)))

# Join on year
al.mld <- winter.wide %>%
  left_join(al.slp %>% select(year, AL.sd), by = "year") %>%
  filter(!is.na(AL.sd))

al.sd.range   <- range(al.mld$AL.sd,   na.rm = TRUE)
goa.ar1.range <- range(al.mld$GOA.ar1, na.rm = TRUE)
ebs.ar1.range <- range(al.mld$EBS.ar1, na.rm = TRUE)

make_mld_dual_plot <- function(dat, ar1.col, ar1.color, region.label, ar1.range) {
  k <- diff(ar1.range) / diff(al.sd.range)
  b <- ar1.range[1] - k * al.sd.range[1]

  dat <- dat %>% filter(!is.na(.data[[ar1.col]]), !is.na(AL.sd))

  r   <- cor(dat[[ar1.col]], dat$AL.sd, use = "complete.obs")
  lbl <- sprintf("r = %.2f", r)

  ggplot(dat, aes(x = year)) +
    geom_line(aes(y = .data[[ar1.col]], color = "MLD AR(1)"), linewidth = 0.7) +
    geom_point(aes(y = .data[[ar1.col]], color = "MLD AR(1)"), size = 1.5) +
    geom_line(aes(y = k * AL.sd + b, color = "AL SLP SD"),
              linewidth = 0.7, linetype = "dashed") +
    geom_point(aes(y = k * AL.sd + b, color = "AL SLP SD"), size = 1.5, shape = 1) +
    annotate("text", x = Inf, y = Inf, label = lbl,
             hjust = 1.1, vjust = 1.5, size = 3.5) +
    scale_y_continuous(
      name     = "MLD AR(1)",
      sec.axis = sec_axis(~ (. - b) / k, name = "AL SLP SD (z-scored)")
    ) +
    scale_color_manual(values = c("MLD AR(1)" = ar1.color, "AL SLP SD" = "gray30")) +
    labs(title = region.label, x = "Year", color = NULL) +
    theme_bw() +
    theme(legend.position   = "bottom",
          axis.title.y.left  = element_text(color = ar1.color),
          axis.text.y.left   = element_text(color = ar1.color),
          axis.title.y.right = element_text(color = "gray30"),
          axis.text.y.right  = element_text(color = "gray30"))
}

p.goa.mld <- make_mld_dual_plot(al.mld, "GOA.ar1", "steelblue4", "GOA", goa.ar1.range)
p.ebs.mld <- make_mld_dual_plot(al.mld, "EBS.ar1", "darkorange4", "EBS", ebs.ar1.range)

p.mld.dual <- p.goa.mld / p.ebs.mld +
  plot_annotation(
    title = "15-yr Rolling AL SLP SD vs. MLD AR(1) — Winter (Nov-Mar)",
    theme = theme(plot.title = element_text(hjust = 0.5))
  )

ggsave("./Figures/AL_SD_MLD_AR1_15yr_dual_axis.png",
       plot = p.mld.dual, width = 10, height = 7, dpi = 150)
message("Saved: Figures/AL_SD_MLD_AR1_15yr_dual_axis.png")

# ---- ERA-SPECIFIC CELLWISE AR(1) MAPS: SST and MLD ----
# Goal: compare AR(1) for era 2 (2003-2024) vs era 1 (1987-2008) at each grid cell.
# Both eras span 22 winters. AR(1) computed from winter (Nov-Mar) annual means.
# Produces a two-panel map: SST AR(1) difference (left), MLD AR(1) difference (right).

era1.yrs <- 1989:2006
era2.yrs <- 2007:2024

# CHECK 1: era lengths
stopifnot(length(era1.yrs) == 18, length(era2.yrs) == 18)
message("Era 1: ", min(era1.yrs), "-", max(era1.yrs),
        " (n = ", length(era1.yrs), " winters)")
message("Era 2: ", min(era2.yrs), "-", max(era2.yrs),
        " (n = ", length(era2.yrs), " winters)")

# Helper: AR(1) from a matrix [cells x time] for a subset of years
compute_era_ar1 <- function(win.mat, all.yrs, era.yrs, min.n = 10) {
  idx <- which(all.yrs %in% era.yrs)
  apply(win.mat[, idx, drop = FALSE], 1, function(v) {
    v <- v[!is.na(v)]
    if (length(v) < min.n) return(NA_real_)
    tryCatch(acf(v, lag.max = 1, plot = FALSE)$acf[2],
             error = function(e) NA_real_)
  })
}

ar1.diff.cache <- "./Output/era_ar1_diff_sst_mld.csv"

if (file.exists(ar1.diff.cache)) {
  message("Loading cached AR(1) difference maps...")
  ar1.diff <- read.csv(ar1.diff.cache)
} else {

  # ---- SST: CDO preprocessing ----
  sst.src        <- "./Data/era5_sst_NP_1950_2025.nc"
  sst.winmon.file <- "./Data/era5_sst_NP_detr_anom_winmon.nc"

  if (!file.exists(sst.winmon.file)) {
    message("CDO: ERA5 SST detrended winter-month anomalies...")
    system(paste(
      cdo, "selmon,11,12,1,2,3 -detrend -ymonsub", sst.src,
      paste0("-ymonmean ", sst.src),
      sst.winmon.file
    ))
  }

  message("Loading ERA5 SST detrended winter anomalies...")
  nc.s     <- nc_open(sst.winmon.file)
  sst.w    <- ncvar_get(nc.s, "sst")        # Kelvin; offset cancels in anomalies
  lons.s   <- ncvar_get(nc.s, "longitude")
  lats.s   <- ncvar_get(nc.s, "latitude")
  # CDO renames valid_time -> time, stored as a dimension variable
  all.names <- c(names(nc.s$var), names(nc.s$dim))
  time.var  <- intersect(c("valid_time", "time"), all.names)[1]
  time.s    <- ncvar_get(nc.s, time.var)
  tun.s     <- ncatt_get(nc.s, time.var, "units")$value
  nc_close(nc.s)

  origin.s    <- sub(".*since ", "", tun.s)
  dates.s     <- as.Date(as.POSIXct(time.s, origin = origin.s, tz = "UTC"))
  years.s     <- as.integer(format(dates.s, "%Y"))
  months.s    <- as.integer(format(dates.s, "%m"))
  winyears.s  <- ifelse(months.s %in% c(11, 12), years.s + 1L, years.s)
  lons.s360   <- ifelse(lons.s < 0, lons.s + 360, lons.s)

  # Crop to NP domain
  lon.idx.s <- which(lons.s360 >= 110 & lons.s360 <= 250)
  lat.idx.s <- which(lats.s   >= 20  & lats.s   <= 66)
  sst.w     <- sst.w[lon.idx.s, lat.idx.s, ]
  lons.s360 <- lons.s360[lon.idx.s]
  lats.s    <- lats.s[lat.idx.s]

  # CHECK 2: SST year coverage spans both eras
  message("SST winter years available: ", min(winyears.s), "-", max(winyears.s))
  stopifnot(all(era1.yrs %in% winyears.s), all(era2.yrs %in% winyears.s))

  nlon.s <- length(lons.s360); nlat.s <- length(lats.s)
  sst.mat <- matrix(sst.w, nrow = nlon.s * nlat.s, ncol = dim(sst.w)[3])
  rm(sst.w); gc()

  all.wyrs.s <- sort(unique(winyears.s))
  message("Computing SST winter means per cell (", nlon.s * nlat.s, " cells)...")
  sst.win <- matrix(NA_real_, nrow = nlon.s * nlat.s, ncol = length(all.wyrs.s))
  for (k in seq_along(all.wyrs.s)) {
    ti <- which(winyears.s == all.wyrs.s[k])
    sst.win[, k] <- rowMeans(sst.mat[, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(sst.mat); gc()

  message("Computing SST AR(1) per cell for each era...")
  sst.ar1.e1   <- compute_era_ar1(sst.win, all.wyrs.s, era1.yrs)
  sst.ar1.e2   <- compute_era_ar1(sst.win, all.wyrs.s, era2.yrs)
  sst.ar1.diff <- sst.ar1.e2 - sst.ar1.e1
  rm(sst.win); gc()

  # CHECK 3: AR(1) values within [-1, 1]
  message("SST AR(1) era1 range: ", paste(round(range(sst.ar1.e1, na.rm=TRUE), 3), collapse=" to "))
  message("SST AR(1) era2 range: ", paste(round(range(sst.ar1.e2, na.rm=TRUE), 3), collapse=" to "))
  message("SST AR(1) diff range: ", paste(round(range(sst.ar1.diff, na.rm=TRUE), 3), collapse=" to "))

  sst.df <- expand.grid(lon = lons.s360, lat = lats.s) %>%
    mutate(ar1.e1   = sst.ar1.e1,
           ar1.e2   = sst.ar1.e2,
           diff     = sst.ar1.diff,
           variable = "SST") %>%
    filter(!is.na(diff))

  # ---- MLD: full NP winmon file ----
  message("Loading ORAS5 MLD detrended winter anomalies (full NP)...")
  nc.m   <- nc_open(mld.winmon.file)
  mld.w  <- ncvar_get(nc.m, "somxl030")
  lons.m <- ncvar_get(nc.m, "nav_lon")
  lats.m <- ncvar_get(nc.m, "nav_lat")
  time.m <- ncvar_get(nc.m, "time_counter")
  tun.m  <- ncatt_get(nc.m, "time_counter", "units")$value
  nc_close(nc.m)

  mld.w[mld.w > 1e10] <- NA
  origin.m   <- sub(".*since ", "", tun.m)
  dates.m    <- as.Date(as.POSIXct(time.m, origin = origin.m, tz = "UTC"))
  years.m    <- as.integer(format(dates.m, "%Y"))
  months.m   <- as.integer(format(dates.m, "%m"))
  winyears.m <- ifelse(months.m %in% c(11, 12), years.m + 1L, years.m)
  lons.m360  <- ifelse(lons.m < 0, lons.m + 360, lons.m)

  # CHECK 4: MLD year coverage
  message("MLD winter years available: ", min(winyears.m), "-", max(winyears.m))
  stopifnot(all(era1.yrs %in% unique(winyears.m)),
            all(era2.yrs[era2.yrs <= max(winyears.m)] %in% unique(winyears.m)))

  nx.m  <- dim(mld.w)[1]; ny.m <- dim(mld.w)[2]; nt.m <- dim(mld.w)[3]
  np.mask.m <- as.vector(lons.m360 >= 110 & lons.m360 <= 250 &
                         lats.m    >= 20  & lats.m    <= 66)
  mld.mat   <- matrix(mld.w, nrow = nx.m * ny.m, ncol = nt.m)
  good.m    <- which(np.mask.m & rowSums(!is.na(mld.mat)) > nt.m * 0.5)
  rm(mld.w); gc()

  all.wyrs.m <- sort(unique(winyears.m))
  message("Computing MLD winter means per cell (", length(good.m), " cells)...")
  mld.win <- matrix(NA_real_, nrow = length(good.m), ncol = length(all.wyrs.m))
  for (k in seq_along(all.wyrs.m)) {
    ti <- which(winyears.m == all.wyrs.m[k])
    mld.win[, k] <- rowMeans(mld.mat[good.m, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mld.mat); gc()

  message("Computing MLD AR(1) per cell for each era...")
  mld.ar1.e1   <- compute_era_ar1(mld.win, all.wyrs.m, era1.yrs)
  mld.ar1.e2   <- compute_era_ar1(mld.win, all.wyrs.m, era2.yrs)
  mld.ar1.diff <- mld.ar1.e2 - mld.ar1.e1

  rm(mld.win); gc()

  # CHECK 5: MLD AR(1) values within [-1, 1]
  message("MLD AR(1) era1 range: ", paste(round(range(mld.ar1.e1, na.rm=TRUE), 3), collapse=" to "))
  message("MLD AR(1) era2 range: ", paste(round(range(mld.ar1.e2, na.rm=TRUE), 3), collapse=" to "))
  message("MLD AR(1) diff range: ", paste(round(range(mld.ar1.diff, na.rm=TRUE), 3), collapse=" to "))

  ar1.e1.full <- rep(NA_real_, nx.m * ny.m)
  ar1.e2.full <- rep(NA_real_, nx.m * ny.m)
  diff.full   <- rep(NA_real_, nx.m * ny.m)
  ar1.e1.full[good.m] <- mld.ar1.e1
  ar1.e2.full[good.m] <- mld.ar1.e2
  diff.full[good.m]   <- mld.ar1.diff

  mld.df <- data.frame(
    lon      = as.vector(lons.m360),
    lat      = as.vector(lats.m),
    ar1.e1   = ar1.e1.full,
    ar1.e2   = ar1.e2.full,
    diff     = diff.full,
    variable = "MLD"
  ) %>%
    filter(!is.na(diff), lon >= 110, lon <= 250, lat >= 20, lat <= 66)

  ar1.diff <- bind_rows(sst.df, mld.df)
  write.csv(ar1.diff, ar1.diff.cache, row.names = FALSE)
  message("Saved: ", ar1.diff.cache)

} # end cache block

# CHECK 6: GOA and EBS area-mean SST AR(1) difference — sign should be positive
#          (era 2 = high AL variability, expect more reddening in GOA/EBS)
goa.check <- ar1.diff %>%
  filter(variable == "SST", lon >= 198, lon <= 231, lat >= 54, lat <= 61)
ebs.check <- ar1.diff %>%
  filter(variable == "SST", lon >= 183, lon <= 203, lat >= 53, lat <= 65)
message("CHECK 6 — SST AR(1) diff in GOA (expect +): ",
        round(mean(goa.check$diff, na.rm = TRUE), 3))
message("CHECK 6 — SST AR(1) diff in EBS (expect +): ",
        round(mean(ebs.check$diff, na.rm = TRUE), 3))

# CHECK 7: spatial coherence — high-coupling MLD cells should mirror SST pattern
mld.goa.check <- ar1.diff %>%
  filter(variable == "MLD", lon >= 198, lon <= 231, lat >= 54, lat <= 61)
message("CHECK 7 — MLD AR(1) diff in GOA (expect consistent sign with SST): ",
        round(mean(mld.goa.check$diff, na.rm = TRUE), 3))

# ---- PLOT ----
# Filter to display domain first so color limits reflect only plotted values
sst.plot.df <- ar1.diff %>% filter(variable == "SST", lon >= 180, !is.na(diff))
mld.plot.df <- ar1.diff %>% filter(variable == "MLD", lon >= 180, !is.na(diff))

col.lim.sst <- max(abs(sst.plot.df$diff), na.rm = TRUE)
col.lim.mld <- max(abs(mld.plot.df$diff), na.rm = TRUE)

make_ar1_diff_map <- function(df, title.label, col.lim) {
  ggplot() +
    geom_point(data = df,
               aes(x = lon, y = lat, color = diff),
               shape = 15, size = 0.3) +
    geom_polygon(data = mapWorld,
                 aes(x = long, y = lat, group = group),
                 fill = "gray80", color = "gray50", linewidth = 0.2) +
    geom_polygon(data = polys,
                 aes(x = lon, y = lat, group = region),
                 fill = NA, color = "black", linewidth = 0.7) +
    scale_color_gradientn(
      colors   = rev(oceColorsPalette(64)),
      limits   = c(-col.lim, col.lim),
      values   = scales::rescale(c(-col.lim, 0, col.lim)),
      na.value = "white",
      name     = "\u0394AR(1)\n(era2\u2212era1)"
    ) +
    coord_cartesian(xlim = c(180, 250), ylim = c(20, 66)) +
    labs(title = title.label, x = "Longitude", y = "Latitude") +
    theme_bw() +
    theme(plot.title      = element_text(hjust = 0.5, size = 9),
          legend.position = "right")
}

p.sst.diff <- make_ar1_diff_map(
  sst.plot.df,
  "SST AR(1): era 2 (2007\u20132024) \u2212 era 1 (1989\u20132006)",
  col.lim.sst
)
p.mld.diff <- make_ar1_diff_map(
  mld.plot.df,
  "MLD AR(1): era 2 (2007\u20132024) \u2212 era 1 (1989\u20132006)",
  col.lim.mld
)

p.ar1.era <- (p.sst.diff | p.mld.diff) +
  plot_annotation(
    title    = "Change in AR(1) between eras: era 2 \u2212 era 1",
    subtitle = paste0("Era 1: 1989\u20132006 (n = 18 winters)  |  ",
                      "Era 2: 2007\u20132024 (n = 18 winters)\n",
                      "Winter = Nov\u2013Mar annual mean; positive = increased reddening in era 2"),
    theme = theme(plot.title    = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 8))
  )

ggsave("./Figures/ERA_AR1_diff_SST_MLD.png",
       plot = p.ar1.era, width = 14, height = 5, dpi = 150)
message("Saved: Figures/ERA_AR1_diff_SST_MLD.png")

# ---- ERA MEAN MLD DIFFERENCE MAP ----
# Same two eras; plot change in mean winter MLD (era2 - era1) instead of AR(1).
# MLD array reloaded independently so this section runs regardless of cache state.

mld.mean.cache <- "./Output/era_mean_mld_diff.csv"

if (file.exists(mld.mean.cache)) {
  message("Loading cached MLD mean difference...")
  mld.mean.df <- read.csv(mld.mean.cache)
} else {
  message("Reloading ORAS5 MLD for era mean computation...")
  nc.mm   <- nc_open(mld.winmon.file)
  mld.mm  <- ncvar_get(nc.mm, "somxl030")
  lons.mm <- ncvar_get(nc.mm, "nav_lon")
  lats.mm <- ncvar_get(nc.mm, "nav_lat")
  time.mm <- ncvar_get(nc.mm, "time_counter")
  tun.mm  <- ncatt_get(nc.mm, "time_counter", "units")$value
  nc_close(nc.mm)

  mld.mm[mld.mm > 1e10] <- NA
  origin.mm   <- sub(".*since ", "", tun.mm)
  dates.mm    <- as.Date(as.POSIXct(time.mm, origin = origin.mm, tz = "UTC"))
  years.mm    <- as.integer(format(dates.mm, "%Y"))
  months.mm   <- as.integer(format(dates.mm, "%m"))
  winyears.mm <- ifelse(months.mm %in% c(11, 12), years.mm + 1L, years.mm)
  lons.mm360  <- ifelse(lons.mm < 0, lons.mm + 360, lons.mm)

  nx.mm <- dim(mld.mm)[1]; ny.mm <- dim(mld.mm)[2]; nt.mm <- dim(mld.mm)[3]
  np.mask.mm <- as.vector(lons.mm360 >= 110 & lons.mm360 <= 250 &
                           lats.mm    >= 20  & lats.mm    <= 66)
  mld.mat.mm <- matrix(mld.mm, nrow = nx.mm * ny.mm, ncol = nt.mm)
  good.mm    <- which(np.mask.mm & rowSums(!is.na(mld.mat.mm)) > nt.mm * 0.5)
  rm(mld.mm); gc()

  all.wyrs.mm <- sort(unique(winyears.mm))
  message("Computing MLD winter means per cell...")
  mld.win.mm <- matrix(NA_real_, nrow = length(good.mm), ncol = length(all.wyrs.mm))
  for (k in seq_along(all.wyrs.mm)) {
    ti <- which(winyears.mm == all.wyrs.mm[k])
    mld.win.mm[, k] <- rowMeans(mld.mat.mm[good.mm, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mld.mat.mm); gc()

  compute_era_mean <- function(win.mat, all.yrs, era.yrs, min.n = 10) {
    idx <- which(all.yrs %in% era.yrs)
    apply(win.mat[, idx, drop = FALSE], 1, function(v) {
      v <- v[!is.na(v)]
      if (length(v) < min.n) return(NA_real_)
      mean(v)
    })
  }

  mld.mean.e1   <- compute_era_mean(mld.win.mm, all.wyrs.mm, era1.yrs)
  mld.mean.e2   <- compute_era_mean(mld.win.mm, all.wyrs.mm, era2.yrs)
  mld.mean.diff <- mld.mean.e2 - mld.mean.e1
  rm(mld.win.mm); gc()

  mean.e1.full  <- rep(NA_real_, nx.mm * ny.mm)
  mean.e2.full  <- rep(NA_real_, nx.mm * ny.mm)
  meandiff.full <- rep(NA_real_, nx.mm * ny.mm)
  mean.e1.full[good.mm]  <- mld.mean.e1
  mean.e2.full[good.mm]  <- mld.mean.e2
  meandiff.full[good.mm] <- mld.mean.diff

  mld.mean.df <- data.frame(
    lon     = as.vector(lons.mm360),
    lat     = as.vector(lats.mm),
    mean.e1 = mean.e1.full,
    mean.e2 = mean.e2.full,
    diff    = meandiff.full
  ) %>%
    filter(!is.na(diff), lon >= 110, lon <= 250, lat >= 20, lat <= 66)

  write.csv(mld.mean.df, mld.mean.cache, row.names = FALSE)
  message("Saved: ", mld.mean.cache)
}

mld.mean.plot.df <- mld.mean.df %>% filter(lon >= 180, !is.na(diff))
col.lim.mean <- max(abs(mld.mean.plot.df$diff), na.rm = TRUE)

p.sst.diff2 <- make_ar1_diff_map(
  sst.plot.df,
  "SST AR(1): era 2 (2007\u20132024) \u2212 era 1 (1989\u20132006)",
  col.lim.sst
)

p.mld.mean.diff <- ggplot() +
  geom_point(data = mld.mean.plot.df,
             aes(x = lon, y = lat, color = diff),
             shape = 15, size = 0.3) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region),
               fill = NA, color = "black", linewidth = 0.7) +
  scale_color_gradientn(
    colors   = rev(oceColorsPalette(64)),
    limits   = c(-col.lim.mean, col.lim.mean),
    values   = scales::rescale(c(-col.lim.mean, 0, col.lim.mean)),
    na.value = "white",
    name     = "\u0394MLD (m)\n(era2\u2212era1)"
  ) +
  coord_cartesian(xlim = c(180, 250), ylim = c(20, 66)) +
  labs(title = "Mean MLD: era 2 (2007\u20132024) \u2212 era 1 (1989\u20132006)",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(plot.title      = element_text(hjust = 0.5, size = 9),
        legend.position = "right")

p.era.mean <- (p.sst.diff2 | p.mld.mean.diff) +
  plot_annotation(
    title    = "Era differences: SST AR(1) change and mean MLD change",
    subtitle = paste0("Era 1: 1989\u20132006 (n = 18 winters)  |  ",
                      "Era 2: 2007\u20132024 (n = 18 winters)\n",
                      "Winter = Nov\u2013Mar annual mean"),
    theme = theme(plot.title    = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 8))
  )

ggsave("./Figures/ERA_SST_AR1_diff_MLD_mean_diff.png",
       plot = p.era.mean, width = 14, height = 5, dpi = 150)
message("Saved: Figures/ERA_SST_AR1_diff_MLD_mean_diff.png")

# ---- SST AR(1) ERA DIFF + SLP PC1-MLD REGRESSION — combined two-panel map ----
# Left:  SST AR(1) era difference (era2 - era1), domain 180-250E, 20-66N
# Right: Cellwise SLP PC1 -> MLD regression over same domain, rerun fresh.
# CDO lon range: -180 to 110 covers 180-250E in -180/180 space.

reg.domain.cache <- "./Output/oras5_mld_slp_pc1_regression_180_250.csv"

if (file.exists(reg.domain.cache)) {
  message("Loading cached 180-250E regression...")
  reg.domain.df <- read.csv(reg.domain.cache)
} else {
  # Reuse the already-cropped NE Pacific file (170-250E, 30-66N) — small enough
  # to load. Mask to lon >= 180 in R. Note: lat coverage is 30-66N.
  mld.ne.file <- "./Data/oras5_mld_NEP_detr_anom_winmon.nc"
  message("Loading MLD anomalies from NE Pacific file, masking to 180-250E...")
  nc.d   <- nc_open(mld.ne.file)
  mld.d  <- ncvar_get(nc.d, "somxl030")
  lons.d <- ncvar_get(nc.d, "nav_lon")
  lats.d <- ncvar_get(nc.d, "nav_lat")
  time.d <- ncvar_get(nc.d, "time_counter")
  tun.d  <- ncatt_get(nc.d, "time_counter", "units")$value
  nc_close(nc.d)

  mld.d[mld.d > 1e10] <- NA
  origin.d   <- sub(".*since ", "", tun.d)
  dates.d    <- as.Date(as.POSIXct(time.d, origin = origin.d, tz = "UTC"))
  years.d    <- as.integer(format(dates.d, "%Y"))
  months.d   <- as.integer(format(dates.d, "%m"))
  winyears.d <- ifelse(months.d %in% c(11, 12), years.d + 1L, years.d)
  lons.d360  <- ifelse(lons.d < 0, lons.d + 360, lons.d)

  domain.mask <- lons.d360 >= 180 & lons.d360 <= 250 & lats.d >= 20 & lats.d <= 66

  pc1.d     <- read.csv("./Output/NP_PC1_winter.csv")
  common.d  <- sort(intersect(unique(winyears.d), pc1.d$year))
  pc1.vec.d <- pc1.d$PC1[match(common.d, pc1.d$year)]

  nx.d <- dim(mld.d)[1]; ny.d <- dim(mld.d)[2]; nt.d <- dim(mld.d)[3]
  mat.d  <- matrix(mld.d, nrow = nx.d * ny.d, ncol = nt.d)
  good.d <- which(domain.mask & rowSums(!is.na(mat.d)) > nt.d * 0.5)
  rm(mld.d); gc()

  message("Computing cellwise winter means (180-250E domain)...")
  win.d <- matrix(NA_real_, nrow = length(good.d), ncol = length(common.d))
  for (k in seq_along(common.d)) {
    ti <- which(winyears.d == common.d[k])
    win.d[, k] <- rowMeans(mat.d[good.d, ti, drop = FALSE], na.rm = TRUE)
  }
  rm(mat.d); gc()

  message("Fitting GLS AR(1) regressions (180-250E domain)...")
  beta.d <- rep(NA_real_, length(good.d))
  pval.d <- rep(NA_real_, length(good.d))
  df.d   <- data.frame(y = NA_real_, x = pc1.vec.d)

  for (k in seq_along(good.d)) {
    y <- win.d[k, ]
    if (sum(!is.na(y)) < 10) next
    df.d$y <- y
    fit <- tryCatch(
      gls(y ~ x, data = df.d,
          correlation = corAR1(form = ~ 1),
          method = "ML", na.action = na.omit),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    s <- summary(fit)$tTable
    if (nrow(s) < 2) next
    beta.d[k] <- s[2, 1]
    pval.d[k] <- s[2, 4]
  }

  beta.full.d <- rep(NA_real_, nx.d * ny.d)
  pval.full.d <- rep(NA_real_, nx.d * ny.d)
  beta.full.d[good.d] <- beta.d
  pval.full.d[good.d] <- pval.d

  reg.domain.df <- data.frame(
    lon  = as.vector(lons.d360),
    lat  = as.vector(lats.d),
    beta = beta.full.d,
    pval = pval.full.d
  ) %>%
    filter(!is.na(beta), lon >= 180, lon <= 250, lat >= 20, lat <= 66)

  write.csv(reg.domain.df, reg.domain.cache, row.names = FALSE)
  message("Saved: ", reg.domain.cache)
}

reg.domain.grid <- reg.domain.df %>%
  mutate(lon = round(lon / 0.5) * 0.5,
         lat = round(lat / 0.5) * 0.5) %>%
  group_by(lon, lat) %>%
  summarise(beta = mean(beta, na.rm = TRUE),
            pval = mean(pval, na.rm = TRUE),
            .groups = "drop") %>%
  complete(lon = seq(180, 250, by = 0.5),
           lat = seq(20,  66,  by = 0.5))

col.lim.reg <- max(abs(reg.domain.df$beta), na.rm = TRUE)

p.reg.domain <- ggplot() +
  geom_point(data = filter(reg.domain.df, !is.na(beta)),
             aes(x = lon, y = lat, color = beta),
             shape = 15, size = 0.5) +
  geom_contour(data = filter(reg.domain.grid, !is.na(pval)),
               aes(x = lon, y = lat, z = pval),
               breaks = 0.05, color = "black", linewidth = 0.4) +
  geom_polygon(data = mapWorld,
               aes(x = long, y = lat, group = group),
               fill = "gray80", color = "gray50", linewidth = 0.2) +
  geom_polygon(data = polys,
               aes(x = lon, y = lat, group = region),
               fill = NA, color = "black", linewidth = 0.7) +
  scale_color_gradientn(
    colors   = rev(oceColorsPalette(64)),
    limits   = c(-col.lim.reg, col.lim.reg),
    values   = scales::rescale(c(-col.lim.reg, 0, col.lim.reg)),
    na.value = "white",
    name     = "\u03b2 (m / PC1)"
  ) +
  coord_cartesian(xlim = c(180, 250), ylim = c(20, 66)) +
  labs(title    = "SLP PC1 \u2192 MLD regression (full period)",
       subtitle = "GLS AR(1); contour = p < 0.05",
       x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(plot.title      = element_text(hjust = 0.5, size = 9),
        plot.subtitle   = element_text(hjust = 0.5, size = 7),
        legend.position = "right")

p.combined <- (p.sst.diff | p.reg.domain) +
  plot_annotation(
    title    = "SST AR(1) era change and SLP PC1\u2192MLD regression",
    subtitle = paste0("Left: SST AR(1) era 2 (2007\u20132024) \u2212 era 1 (1989\u20132006)  |  ",
                      "Right: full-period regression\nDomain: 180\u2013250E, 20\u201366N; ",
                      "winter = Nov\u2013Mar annual mean"),
    theme = theme(plot.title    = element_text(hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5, size = 8))
  )

ggsave("./Figures/ERA_SST_AR1_diff_SLP_MLD_regression.png",
       plot = p.combined, width = 14, height = 5, dpi = 150)
message("Saved: Figures/ERA_SST_AR1_diff_SLP_MLD_regression.png")
